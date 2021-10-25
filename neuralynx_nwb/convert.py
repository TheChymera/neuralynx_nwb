import numpy as np
import pynwb
from datetime import datetime
from dateutil.tz import tzlocal
from os import path
from pynwb import NWBFile
import re
import time

import neo


def reposit_data(
	data_dir='~/.local/share/datalad/',
	data_selection='_vStr_phase_stim/M235/M235-2021-07-16/',
	lab_name='MVDMLab',
	institution='Dartmouth College',
	keywords=[
		'DANDI Pilot',
		],
	experimenter='Manish ...',
	experiment_description= '...',
	debug=True,
	):

	data_dir = path.abspath(path.expanduser(data_dir))
	session_data = path.join(data_dir,data_selection)
	# session_data = '/Users/jimmiegmaz/Desktop/M040-2020-04-28-CDOD11' #for testing

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

	# create a reader
	#reader = neo.io.NeuralynxIO(dirname=session_data) # TODO: newer version should support: , keep_original_times=True)
	reader = neo.io.NeuralynxIO(dirname=session_data, keep_original_times=False) # TODO: newer version should support: , keep_original_times=True)
	reader.parse_header()
	if debug:
		print(reader)

	# seg = reader.read_segment()
	seg = reader.read()

	filename_metadata = re.match(
	'(?P<subject_id>[A-Za-z0-9]*)-(?P<date>20..-..-..)-(?P<task>[A-Za-z]*)(?P<day_of_recording>[0-9]*)$',
	path.basename(session_data)
		).groupdict()
	if debug:
		print(filename_metadata)

	# Those time stamps are in sub-second and not the one we would want to the "session time"
	# time.gmtime(reader.get_event_timestamps()[0][0])
	# TODO: figure out where in this 
	# TODO: figure out what those timestamps in.
	
	

	# Scans through Experimental Keys to extract relevant metadata for NWB file

	## name of ExpKeys file
	keys_name = session_data + '/'  + filename_metadata['subject_id'] + '_' + filename_metadata['date'].replace('-','_') + '_keys.m'

	## read session ExpKeys
	with open (keys_name, 'rt') as keys_file:
		exp_keys = keys_file.read()

	## list of metadata to extract
	metadata_list = ['ExpKeys.species','ExpKeys.hemisphere','ExpKeys.weight','ExpKeys.probeDepth','ExpKeys.target']
	reader = neo.io.NeuralynxIO(dirname=session_data) # TODO: newer version should support: , keep_original_times=True)

	# initialize metadata dictionary
	metadata_keys = dict.fromkeys(metadata_list)

	## extract metadata
	for item in exp_keys.split("\n"):
		for field in metadata_list:
			if field in item:
				metadata_keys[field] = re.search('(?<=\=)(.*?)(?=\;)', item).group(0).strip() 
				metadata_keys[field] = re.sub('[^A-Za-z0-9]+', '', metadata_keys[field])
				if debug:
					print(metadata_keys[field])
				
	## TODO: add surgery details to ExpKeys, including AP and ML coordinates, change probeDepth to mm,
	## add filtering, individual tetrode depth, tetrode referencing

	

	# Metadata which is likely to come from data files and "promotion" metadata records

	# Most likely many could be parsed from the filenames which are likely to encode some of it
	# So "heuristical" converter could establish metadata harvesting from the filenames

	# Session specific
	session_metadata = dict(
		session_id="%(subject_id)s-%(date)s" % filename_metadata,
		session_description="Extracellular ephys recording in the left hemisphere of the nucleus accumbens",  # args[0] in nwbfile
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

	
	# Such NWBFile will be created for each separate file, and then fill up with the corresponding
	#
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

	# add electrode metadata
	# create probe device
	device = nwbfile.create_device(name='silicon probe', description='A4x2-tet-5mm-150-200-121', manufacturer='NeuroNexus')

	# for each channel on the probe
	for chl in reader.header['unit_channels']:
		
		# get tetrode id
		tetrode = re.search('(?<=TT)(.*?)(?=#)', chl[0]).group(0)
		electrode_name = 'tetrode' + tetrode
		
		# get channel id
		channel = re.search('(?<=#)(.*?)(?=#)', chl[0]).group(0)
			   
		if electrode_name not in nwbfile.electrode_groups: # make tetrode if does not exist
		
			description = electrode_name
			location = metadata_keys['ExpKeys.hemisphere'] + ' ' + metadata_keys['ExpKeys.target'] + ' ' + \
				'(' + metadata_keys['ExpKeys.probeDepth'] + ' um)'

			electrode_group = nwbfile.create_electrode_group(electrode_name,
															 description=description,
															 location=location,
															 device=device)
			
		# add channel to tetrode
		nwbfile.add_electrode(id=int(channel),
							x=-1.2, y=float(metadata_keys['ExpKeys.probeDepth']), z=-1.5,
							location=metadata_keys['ExpKeys.target'], filtering='none',
							imp = 0.0, group=nwbfile.electrode_groups[electrode_name])

	

	# append data from different segments

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
	for i, chl in enumerate(reader.header['unit_channels']):
		
		# get tetrode id
		tetrode = re.search('(?<=TT)(.*?)(?=#)', chl[0]).group(0)
		tetrode_name = 'TT' + tetrode
			   
		if tetrode_name not in ephys_waveform.spike_event_series: # make tetrode if does not exist
			
			chl_list = []
			
			for j, group in enumerate(nwbfile.electrodes['group']):
			
				if tetrode in nwbfile.electrodes['group'][j].fields['description']:
					
					chl_list.append(j)
			
			electrode_table_region = nwbfile.create_electrode_table_region(chl_list, tetrode_name)
			
			for s in range(reader.header['nb_segment'][0]):
				
				if s == 0:
					
					waveform = reader.get_spike_raw_waveforms(seg_index=s, unit_index=i)
					
				else:
					
					waveform = np.vstack([waveform,reader.get_spike_raw_waveforms(seg_index=s, unit_index=i)])

			ephys_waveform.create_spike_event_series(tetrode_name,
													 waveform,
													 spk_all[i],
													 electrode_table_region)

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
