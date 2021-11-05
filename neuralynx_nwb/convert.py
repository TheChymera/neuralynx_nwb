import numpy as np
import re
import time
from datetime import datetime
from dateutil.tz import tzlocal
from os import path

import neo
import pynwb
from pynwb import NWBFile
from pynwb.ogen import OptogeneticStimulusSite, OptogeneticSeries
from ndx_optogenetics import OpticFiberImplant, OrthogonalStereotacticTarget


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
	# Experiment specific one
	experiment_metadata = dict(
		experimenter=experimenter,
		experiment_description=experiment_description,
	)

	# create multiple readers, pending resolution of:
	# https://github.com/NeuralEnsemble/python-neo/issues/1042#issuecomment-957297763
	reader_lfp = neo.io.NeuralynxIO(dirname=session_dir,
			keep_original_times=keep_original_times,
			exclude_filename=[
				'WE1.ncs',
				'WE2.ncs',
				'CSC1.ncs',
				'CSC2.ncs',
				'CSC3.ncs',
				'CSC4.ncs',
				'CSC5.ncs',
				'CSC6.ncs',
				'CSC7.ncs',
				'CSC8.ncs',
				'CSC9.ncs',
				'CSC10.ncs',
				'CSC11.ncs',
				'CSC12.ncs',
				'CSC13.ncs',
				'CSC14.ncs',
				'CSC15.ncs',
				'CSC16.ncs',
				'CSC17.ncs',
				'CSC18.ncs',
				'CSC19.ncs',
				'CSC20.ncs',
				'CSC21.ncs',
				'CSC22.ncs',
				'CSC23.ncs',
				'CSC24.ncs',
				'CSC25.ncs',
				'CSC26.ncs',
				'CSC27.ncs',
				'CSC28.ncs',
				'CSC29.ncs',
				'CSC30.ncs',
				'CSC31.ncs',
				'CSC32.ncs',
				],
			)
	reader_csc = neo.io.NeuralynxIO(dirname=session_dir,
			keep_original_times=keep_original_times,
			exclude_filename=[
				'WE1.ncs',
				'WE2.ncs',
				'LFP28.ncs',
				'LFP30.ncs',
				'LFP4.ncs',
				'LFP6.ncs',
				],
			)
	reader_we = neo.io.NeuralynxIO(dirname=session_dir,
			keep_original_times=keep_original_times,
			exclude_filename=[
				'LFP28.ncs',
				'LFP30.ncs',
				'LFP4.ncs',
				'LFP6.ncs',
				'CSC1.ncs',
				'CSC2.ncs',
				'CSC3.ncs',
				'CSC4.ncs',
				'CSC5.ncs',
				'CSC6.ncs',
				'CSC7.ncs',
				'CSC8.ncs',
				'CSC9.ncs',
				'CSC10.ncs',
				'CSC11.ncs',
				'CSC12.ncs',
				'CSC13.ncs',
				'CSC14.ncs',
				'CSC15.ncs',
				'CSC16.ncs',
				'CSC17.ncs',
				'CSC18.ncs',
				'CSC19.ncs',
				'CSC20.ncs',
				'CSC21.ncs',
				'CSC22.ncs',
				'CSC23.ncs',
				'CSC24.ncs',
				'CSC25.ncs',
				'CSC26.ncs',
				'CSC27.ncs',
				'CSC28.ncs',
				'CSC29.ncs',
				'CSC30.ncs',
				'CSC31.ncs',
				'CSC32.ncs',
				],
			)
	reader_lfp.parse_header()
	reader_csc.parse_header()
	reader_we.parse_header()

	reader = reader_csc

	print('Reading from: {}'.format(session_dir))
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
	#     hemisphere=metadata_keys['ExpKeys.hemisphere'],
	#     depth=metadata_keys['ExpKeys.probeDepth'],
	#     region=metadata_keys['ExpKeys.target'],
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
	# create probe device
	device = nwbfile.create_device(name='silicon probe', description='A4x2-tet-5mm-150-200-121', manufacturer='NeuroNexus')

	print(reader.header)

	# for each channel on the probe
	for chl in reader.header['spike_channels']:
		# get tetrode id
		tetrode = re.search('(?<=TT)(.*?)(?=#)', chl[0]).group(0)
		electrode_name = 'tetrode' + tetrode

		# get channel id
		channel = re.search('(?<=#)(.*?)(?=#)', chl[0]).group(0)

		if electrode_name not in nwbfile.electrode_groups: # make tetrode if does not exist
			print('Adding Electrode: {}'.format(electrode_name))
			description = electrode_name
			# # Pending inclusion of ExpKeys data, empty string location for the time being:
			# location_full = metadata_keys['ExpKeys.hemisphere'] + ' ' + metadata_keys['ExpKeys.target'] + ' ' + \
			# 	'(' + metadata_keys['ExpKeys.probeDepth'] + ' um)'
			location_full = ''

			electrode_group = nwbfile.create_electrode_group(electrode_name,
					description=description,
					location=location_full,
					device=device,
					)
			# add channel to tetrode
			# All of these fields should ideally be fetched from ExpKeys fields, and not hard-coded here.
			# Format would be, e.g. `y=float(metadata_keys['ExpKeys.probeDepth'])` or `location=metadata_keys['ExpKeys.target']`
			# This had a higher indent level in the original script, appears that might have been as mistake.
			nwbfile.add_electrode(
					id=int(channel),
					x=-1.2,
					y=2.,
					z=-1.5,
					location='target',
					filtering='none',
					imp = 0.0,
					group=nwbfile.electrode_groups[electrode_name],
					)
	if debug:
		print('Detected the following electrodes: {}'.format(nwbfile.electrode_groups))


	# Start reading actual data, segment-wise
	reader = reader_csc
	seg = reader.read()

	spk_all = []
	wv_all = []
	csc_all_mag = []
	csc_all_time = [];
	beh_all = []

	for s in range(reader.header['nb_segment'][0]):
	#     for i, chl in enumerate(reader.header['unit_channels']):
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


	# add data to nwb file
	from pynwb.ecephys import ElectricalSeries
	from pynwb.ecephys import SpikeEventSeries
	from pynwb.ecephys import EventWaveform
	from pynwb.behavior import BehavioralTimeSeries

	# add .ntt files
	ephys_waveform = EventWaveform()

	# loop through .ntt files
	# print(reader.header)
	for i, chl in enumerate(reader.header['spike_channels']):
		
		# get tetrode id
		tetrode = re.search('(?<=TT)(.*?)(?=#)', chl[0]).group(0)
		tetrode_name = 'TT' + tetrode
		if debug:
			print('Detected tetrode {} from header'.format(tetrode_name))
			   
		if tetrode_name not in ephys_waveform.spike_event_series: # make tetrode if does not exist
			print('Adding Tetrode: {}'.format(tetrode_name))
			chl_list = []
			for j, group in enumerate(nwbfile.electrodes['group']):
				if tetrode in nwbfile.electrodes['group'][j].fields['description']:
					chl_list.append(j)
			
			electrode_table_region = nwbfile.create_electrode_table_region(chl_list, tetrode_name)
			waveform = reader.get_spike_raw_waveforms(spike_channel_index=i)
			'''for s in range(reader.header['nb_segment'][0]):
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
					waveform = np.vstack([waveform,reader.get_spike_raw_waveforms(s, i)])'''
			print(np.shape(waveform))
			print(np.shape(spk_all[i][0]))
			print(spk_all[i])
			print(spk_all[i][0])

			ephys_waveform.create_spike_event_series(tetrode_name,
													 waveform,
													 spk_all[i][0],
													 electrode_table_region,
													 )

	nwbfile.add_acquisition(ephys_waveform)

	# add .ncs files

	chl_list = []

	for chl in reader.header['signal_channels']['id']:
		
		chl_list.append(nwbfile.electrodes['id'][:].index(chl))
		
	electrode_table_region = nwbfile.create_electrode_table_region(chl_list, 'CSC order for time series')

	ephys_ts = ElectricalSeries('CSC data',
								csc_all_mag,
								electrode_table_region,
								timestamps=csc_all_time,
								comments='n/a',
								description='unfiltered CSC data')

	nwbfile.add_acquisition(ephys_ts)

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
