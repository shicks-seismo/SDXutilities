#!/usr/bin/env python
"""
Utilities for working with SDX (Seismic Data Explorer) event files
Functions:
    - sdxtoquakeml: Convert SDX to QuakeML format using ObsPy inventory
      structure

"""
import glob
from itertools import islice
from obspy.core.event import (Catalog, CreationInfo, Event, EventDescription,
                              Origin, ResourceIdentifier, Pick,
                              WaveformStreamID, Arrival, OriginQuality)
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees


def sdxtoquakeml(sdx_dir, out_xml,
                 time_uncertainties=[0.1, 0.2, 0.5, 0.8, 1.5],
                 catalog_description="", catalog_version="",
                 agency_id="", author="", vel_mod_id=""):
    """
    Convert SDX to QuakeML format using ObsPy inventory structure.
    SDX filename prefix is stored under event description.
    Input parameters:
        - sdx_dir: directory containing sdx files (required)
        - out_xml: Filename of quakeML file (required)
        - time_uncertainties: List containing time uncertainities in seconds
          for mapping from weights 0-4, respectively (optional)
        - catalog_description (optional)
        - cat_agency_id (optional)
        - author (optional)
        - vel_mod_id (optional)
    Output:
        - xml catalog in QuakeML format.
    """

    # Prepare catalog
    cat = Catalog(description=catalog_description,
                  creation_info=CreationInfo(
                      author=author, agency_id=agency_id,
                      version=catalog_version))

    # Read in sdx files in directory, recursively
    files = glob.glob("{:}/**/*.sdx".format(sdx_dir), recursive=True)
    if len(files) == 0:
        print("No SDX files found in path. Exiting")
    for sdx_file_path in files:
        print("Working on ", sdx_file_path.split('/')[-1])

        # Set-up event
        evt_id = (sdx_file_path.split('/')[-1])[:-4]
        event = Event(event_type="earthquake", creation_info=CreationInfo(
            author=author, agency_id=agency_id),
            event_descriptions=[EventDescription(text=evt_id)])

        # Get station details, append to arrays
        sdx_file = open(sdx_file_path, "r")
        stations = []
        for line in sdx_file:
            if line.rstrip() == "station":
                sdxstation = list(islice(sdx_file, 5))
                stations.append([sdxstation[1].split()[0],
                                 float(sdxstation[2].split()[0]),
                                 float(sdxstation[3].split()[0]),
                                 float(sdxstation[4].split()[0])])
        sdx_file.close()

        # Find origin details, append to origin object
        sdx_file = open(sdx_file_path, "r")
        found_origin = False
        for line in sdx_file:
            if line.rstrip() == "origin":
                found_origin = True
                sdxorigin = list(islice(sdx_file, 17))
                orig_time = ("{:}T{:}".
                             format(sdxorigin[1][0:10].replace(".", "-"),
                                    sdxorigin[1][11:23]))
                evt_lat = float(sdxorigin[2].split()[0])
                evt_lon = float(sdxorigin[3].split()[0])
                evt_depth = float(sdxorigin[4].split()[0])
                creation_time = UTCDateTime(
                    "{:}T{:}".format(sdxorigin[16].split()[6][0:10]
                                     .replace(".", "-"),
                                     sdxorigin[16].split()[6][11:23]))
                num_arrivals = int(sdxorigin[12].split()[0])
                num_arrivals_p = (int(sdxorigin[12].split()[0]) -
                                  int(sdxorigin[12].split()[1]))
                min_dist = float(sdxorigin[12].split()[9])
                max_dist = float(sdxorigin[12].split()[10])
                med_dist = float(sdxorigin[12].split()[11])
                max_az_gap = float(sdxorigin[12].split()[6])

                origin = Origin(time=UTCDateTime(orig_time), longitude=evt_lon,
                                latitude=evt_lat, depth=evt_depth*-1000,
                                earth_model_id=vel_mod_id,
                                origin_type="hypocenter",
                                evaluation_mode="manual",
                                evaluation_status="confirmed",
                                method_id=ResourceIdentifier(id="SDX_hypo71"),
                                creation_info=CreationInfo(
                                    creation_time=creation_time, author=author,
                                    agency_id=agency_id),
                                quality=OriginQuality(
                                    associated_phase_count=num_arrivals,
                                    used_phase_count=num_arrivals,
                                    associated_station_count=num_arrivals_p,
                                    used_station_count=num_arrivals_p,
                                    azimuthal_gap=max_az_gap,
                                    minimum_distance=min_dist,
                                    maximum_distance=max_dist,
                                    median_distance=med_dist))
                event.origins.append(origin)
        sdx_file.close()

        # Skip event if no computed origin
        if found_origin is False:
            print("No origin found ... skipping event")
            continue

        # Get pick details, append to pick and arrival objects
        sdx_file = open(sdx_file_path, "r")
        found_pick = False
        for line in sdx_file:
            if line.rstrip() == "pick":
                found_pick = True
                sdxpick = list(islice(sdx_file, 15))
                pick_time = UTCDateTime(
                    "{:}T{:}".format(sdxpick[1][0:10].replace(".", "-"),
                                     sdxpick[1][11:23]))
                network = sdxpick[2].split()[0]
                station = sdxpick[2].split()[1]
                location = sdxpick[2].split()[2]
                if "NOT_SET" in location:
                    location = ""
                channel = sdxpick[2].split()[3]
                onset = sdxpick[8].split()[0]
                if onset == "0":
                    pickonset = "emergent"
                elif onset == "1":
                    pickonset = "impulsive"
                elif onset == "2":
                    pickonset = "questionable"
                phase = sdxpick[9].split()[0]
                polarity = sdxpick[10].split()[0]
                if polarity == "0":
                    pol = "positive"
                elif polarity == "1":
                    pol = "negative"
                elif polarity == "2":
                    pol = "undecidable"
                weight = int(sdxpick[11].split()[0])
                creation_time = UTCDateTime(
                    "{:}T{:}".format(sdxpick[14].split()[6][0:10]
                                     .replace(".", "-"),
                                     sdxpick[14].split()[6][11:23]))
                pick = Pick(time=pick_time,
                            waveform_id=WaveformStreamID(
                                network_code=network, station_code=station,
                                location_code=location, channel_code=channel),
                            time_errors=time_uncertainties[weight],
                            evaluation_mode="manual",
                            evaluation_status="confirmed", onset=pickonset,
                            phase_hint=phase, polarity=pol,
                            method_id=ResourceIdentifier(id="SDX"),
                            creation_info=CreationInfo(
                                creation_time=creation_time))
                event.picks.append(pick)

                # Compute azimuth, distance, append to arrival object
                for i in range(0, len(stations)):
                    if stations[i][0] == station:
                        azimuth = (gps2dist_azimuth(evt_lat, evt_lon,
                                                    stations[i][1],
                                                    stations[i][2])[1])
                        dist_deg = locations2degrees(evt_lat, evt_lon,
                                                     stations[i][1],
                                                     stations[i][2])
                        arrival = Arrival(phase=phase,
                                          pick_id=pick.resource_id,
                                          azimuth=azimuth, distance=dist_deg,
                                          time_weight=1.00)
                        event.origins[0].arrivals.append(arrival)

        # Skip event if no picks
        if found_pick is False:
            print("No picks found ... skipping event")
            continue

        # Append event to catalogue
        event.preferred_origin_id = sdx_file_path.split("/")[-1][:-4]
        cat.events.append(event)

        sdx_file.close()

    cat.write(out_xml, format="QUAKEML")
