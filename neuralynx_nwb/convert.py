import numpy as np
import re
import time
from copy import deepcopy
from datetime import datetime
from dateutil.tz import tzlocal
from os import path, listdir

import pynwb
from pynwb import NWBFile
from pynwb.ogen import OptogeneticStimulusSite, OptogeneticSeries
from ndx_optogenetics import OpticFiberImplant, OrthogonalStereotacticTarget


def _create_neuralynx_group_readers(session_dir, debug=False, keep_original_times=False):
	"""
	Create reader objects for different section structures.

	Parameters
	----------
	sessions_dir : str
		Session directory for the data in which to construct readers.
	debug : bool, optional
		Whether to print debugging output to the console.
	keep_original_times : bool, optional
		Whether to keep the original times of the data arrays, this argument is passed to `neo.io.NeuralynxIO()`

	Returns
	-------
	readers : dict containing keys which are strings and values which are `neo.io.NeuralynxIO` objects
		Thie `neo.io.NeuralynxIO` objects in this class support NCS, NEV, NSE and NTT file formats.
			* NCS contains signals for one channel
			* NEV contains events
			* NSE contains spikes and waveforms for mono electrodes
			* NTT contains spikes and waveforms for tetrodes
	
	Notes
	-----
	This function may become deprecated pending resolution of:
		https://github.com/NeuralEnsemble/python-neo/issues/1042#issuecomment-957297763
	"""

	import neo
	
	print('Reading from: {}'.format(session_dir))
	
	files_dict = {}
	for i_file in listdir(session_dir):
		try:
			file_prefix = re.findall('(?P<prefix>[A-Za-z]*)[0-9]*\.ncs$', i_file)[0]
		except IndexError:
			pass
		else:
			if not file_prefix in files_dict.keys():
				files_dict[file_prefix] = []
			files_dict[file_prefix].append(i_file)
	if debug:
		print('Created the following dictionary of files based on the {} session directory:'.format(session_dir))
		print(files_dict)

	readers = {}
	for prefix in files_dict.keys():
		exclude_dict = deepcopy(files_dict)
		exclude_dict.pop(prefix)
		exclude_list = exclude_dict.values()
		exclude_list = [val for sublist in exclude_list for val in sublist]
		reader = neo.io.NeuralynxIO(dirname=session_dir,
				keep_original_times=keep_original_times,
				exclude_filename=exclude_list,
				)
		reader.parse_header()
		readers[prefix] = reader
		if debug:
			textfile = open('{}_reader_header.log'.format(prefix), "w")
			textfile.write('{}\n'.format(str(reader.header)))
			textfile.close()

	return readers

def _read_data_segments(reader, debug=False):
	"""
	Return NumPy arrays with data segments from a `neo.io.NeuraLynxIO` reader object.

	Parameters
	----------
	reader : neo.io.NeuralynxIO
		Reader object for which to read data segments
	debug : bool, optional
		Whether to print debugging output to the console.

	Returns
	-------
	spk_all : numpy.array 
		read in spike data
	wv_all
	csc_all_mag
	csc_all_times 
	beh_all
	"""

	# reader at this point needs to be CSC, e.g. `reader = readers['CSC']`
	seg = reader.read()

	spk_all = []
	wv_all = []
	csc_all_mag = []
	csc_all_time = [];
	beh_all = []

	# wv and spk might need to be parsed from another reader (check shapes printed at the end).
	for s in range(reader.header['nb_segment'][0]):
	#	 for i, chl in enumerate(reader.header['unit_channels']):
		if s == 0:
			for i, chl in enumerate(seg[0].segments[s].spiketrains):
				spk_all.append([seg[0].segments[s].spiketrains[i].times])
				
			csc_all_mag = seg[0].segments[s].analogsignals[0].magnitude
			csc_all_time = seg[0].segments[s].analogsignals[0].times
			
			for i, chl in enumerate(seg[0].segments[s].events):
				beh_all.append([seg[0].segments[s].events[i].times])
				
		else:
			for i, chl in enumerate(seg[0].segments[s].spiketrains):
				spk_all[i] = np.append(spk_all[i],[seg[0].segments[s].spiketrains[i].times])
				
			csc_all_mag = np.vstack([csc_all_mag,seg[0].segments[s].analogsignals[0].magnitude])
			csc_all_time = np.append(csc_all_time,seg[0].segments[s].analogsignals[0].times)
			
			for i, chl in enumerate(seg[0].segments[s].events):
				beh_all[i] = np.append(beh_all[i],seg[0].segments[s].events[i].times)

	if debug:
		print(f"spk_all shape: {np.shape(spk_all)}")
		print(f"wv_all shape: {np.shape(wv_all)}")
		print(f"csc_all_mag shape: {np.shape(csc_all_mag)}")
		print(f"csc_all_time shape: {np.shape(csc_all_time)}")
		print(f"beh_all shape: {np.shape(beh_all)}")
	
	return spk_all, wv_all, csc_all_mag, csc_all_time, beh_all

def _setup_electrodes(reader, nwbfile, device_name=None, debug=False):
	"""
	Add electrode groups and continuous signal channels to `pynwb.NWBFile` object while determining the correct assignment.

	Parameters
	----------
	reader : neo.io.NeuralynxIO
		Reader object covering continuous signal channels ("ncs" extension always, "CSC" prefix as per Manish convention) and spike channels ("ntt" extension always, "TT" prefix as per Manish convention).
	nwbfile : pynwb.NWBFile
		NWB file to which to add metadata.
	device_name : str, optional
		Name of the device which to assign to the channels.
		This parameter is not considered if the nwbfile passed only contains one device.
	debug : bool, optional
		Whether to print debugging output to the console.

	Returns
	-------
	pynwb.NWBfile
		NWB file object with added electrodes and electrode groups.
	"""

	# Assign device
	if len(nwbfile.devices.keys()) == 1:
		device = list(nwbfile.devices.values())[0]
	elif isinstance(device_name, str):
		device = nwbfile.devices[device_name]
	else:
		raise ValueError('Please specify a string `device name` parameter to `_setup_electrodes()`.')

	# Set up channels
	electrode_groups = {}
	for chl in reader.header['spike_channels']:
		electrode_matching = '^ch(?P<electrode_group>[a-z,A-Z,0-9]*?)#(?P<channel_nr>[0-9]*?)#.*?$'
		channel_info = re.search(electrode_matching, chl['name']).groupdict()
		electrode_group_name = channel_info['electrode_group']
		channel_nr = int(channel_info['channel_nr'])
		if channel_nr not in electrode_groups.keys():
			electrode_groups[channel_nr] = electrode_group_name
		if electrode_group_name not in nwbfile.electrode_groups: # make tetrode if does not exist
			if debug:
				print('Adding Electrode Group: {}'.format(electrode_group_name))
			description = electrode_group_name
			# # Pending inclusion of ExpKeys data, empty string location for the time being:
			# location_full = metadata_keys['ExpKeys.hemisphere'] + ' ' + metadata_keys['ExpKeys.target'] + ' ' + \
			# 	'(' + metadata_keys['ExpKeys.probeDepth'] + ' um)'
			location_full = ''
			electrode_group = nwbfile.create_electrode_group(electrode_group_name,
					description=description,
					location=location_full,
					device=device,
					)
	if debug:
		print('\nDetected the following electrode groups:\n{}'.format(nwbfile.electrode_groups))

	for channel_nr in reader.header['signal_channels']['id']:
		channel_nr = int(channel_nr)
		electrode_group_name = electrode_groups[channel_nr]
		# add channel to tetrode
		# All of these fields should ideally be fetched from ExpKeys fields, and not hard-coded here.
		# Format would be, e.g. `y=float(metadata_keys['ExpKeys.probeDepth'])` or `location=metadata_keys['ExpKeys.target']`
		# This had a higher indent level in the original script, appears that might have been as mistake.
		if debug:
			print('Adding Signal Channel: {}'.format(channel_nr))
		nwbfile.add_electrode(
			id=channel_nr,
			x=-1.2,
			y=2.,
			z=-1.5,
			location='target',
			filtering='none',
			imp = 0.0,
			group=nwbfile.electrode_groups[electrode_group_name],
			)

	return nwbfile


def reposit_data(
	data_dir='~/.local/share/datalad/',
	data_selection='vStr_phase_stim/M235/M235-2021-07-16/',
	lab_name='MVDMLab',
	institution='Dartmouth College',
	keywords=[
		'DANDI Pilot',
		],
	experimenter='Manish Mohapatra',
	experiment_description='...',
	debug=True,
	session_description='Extracellular ephys recording in the ventral Striatum',
	keep_original_times=True,
	):

	data_dir = path.abspath(path.expanduser(data_dir))
	session_dir = path.normpath(path.join(data_dir,data_selection))

	# Common lab wide metadata
	lab_metadata = dict(
		lab=lab_name,
		institution=institution,
		keywords=keywords,
	)
	# Experiment specific metadata
	experiment_metadata = dict(
		experimenter=experimenter,
		experiment_description=experiment_description,
	)

	readers = _create_neuralynx_group_readers(session_dir, debug=debug, keep_original_times=keep_original_times)

	filename_metadata = re.match(
		'(?P<subject_id>[A-Za-z0-9]*)-(?P<date>20..-..-..)$',
		path.basename(session_dir),
		).groupdict()
	if debug:
		print('Acquired the following metadata from path: {}'.format(filename_metadata))

	# Those time stamps are in sub-second and not the one we would want to the "session time"
	# time.gmtime(reader.get_event_timestamps()[0][0])
	# TODO: figure out where in this 
	# TODO: figure out what those timestamps in.


	# Scans through Experimental Keys to extract relevant metadata for NWB file

	## name of ExpKeys file
	keys_filename = filename_metadata['subject_id'] + '_' + filename_metadata['date'].replace('-','_') + '_keys.m'
	keys_path = path.join(session_dir,keys_filename)

	#  ## read session ExpKeys
	#  with open (keys_path, 'rt') as keys_file:
	#  	exp_keys = keys_file.read()

	## list of metadata to extract and initialize dictionary
	metadata_list = ['ExpKeys.species','ExpKeys.hemisphere','ExpKeys.weight','ExpKeys.probeDepth','ExpKeys.target']
	metadata_keys = dict.fromkeys(metadata_list)

	#  ## extract metadata
	#  for item in exp_keys.split("\n"):
	#  	for field in metadata_list:
	#  		if field in item:
	#  			metadata_keys[field] = re.search('(?<=\=)(.*?)(?=\;)', item).group(0).strip() 
	#  			metadata_keys[field] = re.sub('[^A-Za-z0-9]+', '', metadata_keys[field])
	#  			if debug:
	#  				print(metadata_keys[field])
	#  			
	#  ## TODO: add surgery details to ExpKeys, including AP and ML coordinates, change probeDepth to mm,
	#  ## add filtering, individual tetrode depth, tetrode referencing
	
	# Metadata which is likely to come from data files and "promotion" metadata records
	# Most likely many could be parsed from the filenames which are likely to encode some of it
	# So "heuristical" converter could establish metadata harvesting from the filenames

	# Session specific
	session_metadata = dict(
		session_id="%(subject_id)s-%(date)s" % filename_metadata,
		session_description=session_description,
		session_start_time=datetime.now(tzlocal()), # TEMP  # args[2] in nwbfile; TODO needs to be datetime
	)
	subject_metadata = dict(
		subject_id=filename_metadata['subject_id'],
		weight=metadata_keys['ExpKeys.weight'],
		age="TODO",  # duplicate with session_start_time and date_of_birth but why not?
		species=metadata_keys['ExpKeys.species'],
		sex="female",
	#	 hemisphere=metadata_keys['ExpKeys.hemisphere'],
	#	 depth=metadata_keys['ExpKeys.probeDepth'],
	#	 region=metadata_keys['ExpKeys.target'],
		date_of_birth=datetime.now(tzlocal()), # TEMP: TODO
	)
	surgery_metadata = dict(
		surgery="Headbar on xx/xx/2020, craniotomy over right hemisphere on xx/xx/2020, craniotomy over left hemisphere on xx/xx/2020. All surgeries performed by JG."
	)
	# Actually probably only "identifier" should be file specific, the rest common across files
	# we would like to produce: separate for .ncs, .ntt, behavioral metadata, etc
	file_metadata = dict(
		source_script="somescript-not-clear-whyneeds to be not empty if file_name is provided",
		source_script_file_name="TODO", # __file__,
	)

	# common filename prefix - let's mimic DANDI filenaming convention right away
	filename_prefix = "sub-{subject_id}_ses-{session_id}".format(**subject_metadata, **session_metadata)
	# the rest will be specific to the corresponding file. E.g. we will have separate
	#  - `_probe-<name>_ecephys.nwb` (from each .ncs) - contineous data from each tetrode. probably chunked and compressed
	#  - `_???_ecephys.nwb` (from each .ntt) - spike detected windowed data. 
	#  - `_behav.mpg` + `_behav.nwb` - video recording and metadata (including those .png?) for behavior component within experiment recording session
	# Pretty much we need to establish a framework where EVERY file present would be
	# provided

	if debug:
		print(subject_metadata)

	# TODO add event labels to metadata
	
	# Such NWBFile will be created for each separate file, and then fill up with the corresponding data
	filename_suffix = "TODO"
	nwbfile = NWBFile(
		identifier="{}_{}".format(filename_prefix, filename_suffix), # args[1] in nwbfile, may be just UUID? not sure why user has to provide it really
		subject=pynwb.file.Subject(**subject_metadata),
		**lab_metadata,
		**experiment_metadata,
		**session_metadata,
		**surgery_metadata,
		**file_metadata,)

	if debug:
		print(nwbfile.identifier)

	# Add electrode metadata
	# Create probe device
	nwbfile.create_device(name='silicon probe', description='A4x2-tet-5mm-150-200-121', manufacturer='NeuroNexus')

	reader = readers['CSC']
	nwbfile = _setup_electrodes(reader, nwbfile, debug=debug)
	# This doesn't list the signal channels for some reason
	#if debug:
	#	print('Detected the following channels: {}'.format(nwbfile.electrodes))

	spk_all, wv_all, csc_all_mag, csc_all_time, beh_all = _read_data_segments(reader, debug=debug)

	# add data to nwb file
	from pynwb.ecephys import ElectricalSeries
	from pynwb.ecephys import SpikeEventSeries
	from pynwb.ecephys import EventWaveform
	from pynwb.behavior import BehavioralTimeSeries

	# Add continuous signal data from .ncs files;
	# reader needs to be CSC for this.
	chl_list = []

	for chl in reader.header['signal_channels']:
		# reference by identifier and not by numbered index in the electrode listing as well
		channel_nr = int(chl['id'])
		chl_list.append(nwbfile.electrodes['id'][:].index(channel_nr))
		
	electrode_table_region = nwbfile.create_electrode_table_region(chl_list, 'CSC order for time series')

	ephys_ts = ElectricalSeries('CSC data',
		csc_all_mag,
		electrode_table_region,
		timestamps=csc_all_time,
		comments='n/a',
		description='unfiltered CSC data',
		)

	nwbfile.add_acquisition(ephys_ts)


	# Add spike data from .ntt files
	ephys_waveform = EventWaveform()

	# print(reader.header)
	for i, chl in enumerate(reader.header['spike_channels']):
		electrode_matching = '^ch(?P<electrode_group>[a-z,A-Z,0-9]*?)#(?P<channel_nr>[0-9]*?)#.*?$'
		channel_info = re.search(electrode_matching, chl['name']).groupdict()
		electrode_group_name = channel_info['electrode_group']
		channel_nr = int(channel_info['channel_nr'])
		if debug:
			print('Detected tetrode {} from header'.format(electrode_group_name))

		# It should be in since we just created them based on the selfsame strings.
		#if electrode_group_nr not in ephys_waveform.spike_event_series: # make tetrode if does not exist
		#	print('Adding Tetrode: {}'.format(tetrode_name))
		#	chl_list = []
		#	for j, group in enumerate(nwbfile.electrodes['group']):
		#		if tetrode in nwbfile.electrodes['group'][j].fields['description']:
		#			chl_list.append(j)
			
		#electrode_table_region = nwbfile.create_electrode_table_region(chl_list, tetrode_name)
		waveform = reader.get_spike_raw_waveforms(spike_channel_index=i)
		print("waveform shape", np.shape(waveform))
		print("reader header nb_segment", reader.header['nb_segment'])
		print("reader header nb_segment [0]", reader.header['nb_segment'][0])
		for s in range(reader.header['nb_segment'][0]):
			print('s = {}, i = {}'.format(s,i))
			if s == 0:
				# Pending: https://github.com/NeuralEnsemble/python-neo/issues/1046
				#waveform = reader.get_spike_raw_waveforms(seg_index=s, unit_index=i)
				#waveform = reader.get_spike_raw_waveforms(seg_index=s)
				waveform = reader.get_spike_raw_waveforms(s, i)
			else:
				# Pending: https://github.com/NeuralEnsemble/python-neo/issues/1046
				#waveform = np.vstack([waveform,reader.get_spike_raw_waveforms(seg_index=s, unit_index=i)])
				#waveform = reader.get_spike_raw_waveforms(seg_index=s)
				waveform = np.vstack([waveform,reader.get_spike_raw_waveforms(s, i)])
		print(np.shape(spk_all[i][0]))
		print(spk_all[i])
		print(spk_all[i][0])

		ephys_waveform.create_spike_event_series(electrode_group_name,
			waveform,
			spk_all[i][0],
			electrode_table_region,
			)

	nwbfile.add_acquisition(ephys_waveform)


	# nwbfile.add_unit(id=1, electrodes=[0])
	# nwbfile.add_unit(id=2, electrodes=[0])

	# add .evt file

	beh_ts = BehavioralTimeSeries()

	# loop through events
	for i, chl in enumerate(reader.header['event_channels']):
			
		beh_ts.create_timeseries(str(chl),
									timestamps = beh_all[i])

	nwbfile.add_acquisition(beh_ts)
	
	# Save the generated file
	from pynwb import NWBHDF5IO

	# TODO: I think we should right away use dandi-cli provided API to create the filename based on metadata
	# in the NWBFile
	with NWBHDF5IO('BCD_example.nwb', 'w') as io:
		io.write(nwbfile)

	# Optogenetic Stimulation:
	ogs_site = OptogeneticStimulusSite(
		name='TODO',
		device='TODO',
		description='TODO',
		excitation_lambda='TODO',
		location='TODO',
		)
	ogs_series = OptogeneticSeries(
		name='TODO',
		data='TODO',
		site=ogs_site,
		resolution='TODO',
		conversion='TODO',
		timestampe='TODO',
		starting_time='TODO',
		rate='TODO',
		comments='TODO',
		description='TODO',
		control='TODO',
		control_description='TODO',
		)
