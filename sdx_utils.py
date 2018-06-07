#!/usr/bin/env python
"""
Utilities for working with SDX (Seismic Data Explorer) event files
Functions:
    - sdxtoquakeml: Convert SDX to QuakeML format using ObsPy inventory
      structure.
    - quakemltosdx: Convert QuakeML catalogue into event SDX files using ObsPy
      inventory structure.

"""
import glob
from itertools import islice
import numpy as np
from obspy.core.event import (Catalog, CreationInfo, Event, EventDescription,
                              Origin, ResourceIdentifier, Pick,
                              WaveformStreamID, Arrival, OriginQuality)
from obspy import UTCDateTime, read_events
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

        # Set preferred origin and append event to catalogue
        event.preferred_origin_id = event.origins[0].resource_id
        cat.events.append(event)

        sdx_file.close()

    cat.write(out_xml, format="QUAKEML")


def quakemltosdx(xml_in, station_file, loc_param=[],
                 mseed_list_type=0, mseed_path="",
                 time_uncertainties=[0.1, 0.2, 0.5, 0.8, 1.5]):
    """
    Python sript to convert quakeml file to SDX (Seismic Data Explorer) format
    Usage: QUAKEML_to_SDX.py "[xml_filename] [station_file]"
    Input parameters:
        - xml_in: File containing catalogue in QuakeML format
        - station_file: Text file containing station locations (station name,
            latitude (dd), longitude (dd), elevation (m) - space delimited
        - loc_param: List containing relocation parameters:
            [Path and filename of 1-D velocity model, Vp/Vs_ratio,
             depth_start, depth_step, num_steps, x_near, x_far]
            (optional)
            - mseed_list_type: Trace listing type:
                -0 = One multiplexed miniSEED for each event (default) or
                -1 = one miniSEED file per channel per station.
                (not currently implemented)
        - mseed_path: Path containing miniSEED files for each trace (optional)
        - time_uncertainties: List containing time uncertainities in seconds
            for mapping from weights 0-4, respectively (optional)
            TODO: Choose whether single mSEED file or full directory.
    """

    # Read in catalogue, loop over each event,
    # and define event name based on origin time
    cat = read_events(xml_in)
    for n_evt, event in enumerate(cat):

        # Append origin details if existing
        if event.preferred_origin() is not None:
            origin = event.preferred_origin()
        elif event.preferred_origin() is None and len(event.origins) > 0:
            origin = event.origins[0]

        # Get event ID from description. If no description get from origin
        # time. If no origin, get from first pick time.
        if (event.event_descriptions[0].text is not None and
                event.event_descriptions[0].text[0] == "e"):
            ev_id = event.event_descriptions[0].text
        elif (event.event_descriptions[0].text is None or
              event.event_descriptions[0].text[0] != "e"
              and origin is not None):
            ev_id = ("e{:}{:02g}{:02g}.{:02g}{:02g}{:02.0f}".format(
                origin.time.year, origin.time.month, origin.time.day,
                origin.time.hour, origin.time.minute, origin.time.second))
        elif (event.event_descriptions[0].text is None and origin is None and
              len(event.picks) > 0):
            pick = event.picks[0]
            ev_id = ("e{:}{:02g}{:02g}.{:02g}{:02g}{:02.0f}".format(
                pick.time.year, pick.time.month, pick.time.day,
                pick.time.hour, pick.time.minute, pick.time.second))
        else:
            print("Cannot determine an event ID, skipping event ", n_evt)
            continue

        print("Working on event: ", ev_id)

        # Open SDX file for writing
        out_sdx = "{:}.sdx".format(ev_id)
        w = open(out_sdx, "w")

        # Write header info
        w.write("[+VERSIONSET]\n2\n0\n0\nbasil\n[-VERSIONSET]\n")
        w.write("[+INFOSET]\n\n[-INFOSET]\n")

        # Write synthetics / relocation info (if given)
        w.write("[+SYNTHETICSSET]\n")
        if loc_param:
            w.write("{:}\n{:}\n{:}\n{:}\n{:}\n{:}\n{:}\nF\n\n".format(
                loc_param[1], loc_param[0], loc_param[2], loc_param[3],
                loc_param[4], loc_param[5], loc_param[6]))
        w.write("[-SYNTHETICSSET]\n")

        # Write trace files (if given)
        # TODO: Allow for multiple msd files in single directory.
        w.write("[+TRACESET]\n")
        if mseed_path:
            w.write("{:}\nMS|{:}.msd\n".format(mseed_path, ev_id))
        w.write("[-TRACESET]\n")
        w.write("[+ARCLINKTRACESET]\n[-ARCLINKTRACESET]\n")
        w.write("[+DBTRACESET]\n[-DBTRACESET]\n")

        # Now insert picks
        w.write("[+PICKSET]\n")
        picked_stations = []

        # Loop over arrivals and find corresponding picks
        # TODO: Include non-associated picks?
        for arrival in origin.arrivals:
            w.write("pick\n")
            for _pick in event.picks:
                if arrival.pick_id == _pick.resource_id:
                    pick = _pick

            # Define SNCL info
            network = pick.waveform_id.network_code
            station = pick.waveform_id.station_code
            channel = pick.waveform_id.channel_code

            # Put any hydrophone picks back to the vertical to ensure no
            # ambiguity in channel names and what SDX can read in.
            # This is a PILAB-specific requirement.
            if channel == "BDH":
                channel = "BHZ"
            elif channel == "EDH":
                channel = "HHZ"

            location = pick.waveform_id.location_code
            if location == "" or location is None or location == "None":
                location_full = "<NOT_SET>"
                location_short = ""
            else:
                location_full = location
                location_short = location

            # Define pick time, pick creation time and pick ID, then add to
            # list of picks associated to event/origin
            pick_time = (str(pick.time.datetime).replace(":", "")
                         .replace("-", "").replace(".", "").replace(" ", "")
                         [:-1])
            pick_create_time = (str(pick.creation_info.creation_time).
                                replace("-", ".").replace(" ", "-").
                                replace("T", "-")[:-4])
            pick_id = ("{:}-{:}.{:}.{:}.{:}".format(
                pick_time, network, station, location_short, channel))
            picked_stations.append(pick_id)
            pick_time = (str(pick.time.datetime)
                         .replace("-", ".").replace(" ", "-")[:-3])

            # Define picking uncertainty
            if pick.time_errors["uncertainty"] is not None:
                pick_time_error = pick.time_errors["uncertainty"]
                # Find pick weighting (0-4) corresponding to uncertainity based
                # on index of time_uncertainties list
                pick_weight = (np.where(
                    np.array(time_uncertainties) >= pick_time_error)[0][0])
            else:
                # Set defaults if not existing
                pick_time_error = "0"
                pick_weight = "0"

            # Define ray geometry values
            if pick.backazimuth is None:
                baz = "0"  # Default
            else:
                baz = pick.backazimuth
            if pick.horizontal_slowness is None:
                hsl = "0"  # Default
            else:
                hsl = pick.horizontal_slowness

            # Define onset type and polarity
            if pick.onset == "emergent":
                pick_onset = "0"
            elif pick.onset == "impulsive":
                pick_onset = "1"
            elif pick.onset == "questionable":
                pick_onset = "2"
            else:
                pick_onset = "2"  # Default
            if pick.polarity == "positive":
                pick_polarity = "0"
            elif pick.polarity == "negative":
                pick_polarity = "1"
            elif pick.polarity == "undecidable":
                pick_polarity = "2"
            else:
                pick_polarity = "2"  # Default

            # Now write the pick details
            w.write("{:}\n".format(pick_id))
            w.write("{:} {:} -99999 -99999 -99999\n".format(
                pick_time, pick_time_error))
            w.write("{:} {:} {:} {:} <NOT_SET>\n".format(
                network, station, location_full, channel))
            w.write("1-20%2C%202%2C%20Z\n<NOT_SET>\n"
                    "{:} -99999 -99999 -99999 -99999\n"
                    "{:} -99999 -99999 -99999 -99999\n<NOT_SET>\n"
                    .format(baz, hsl))
            w.write("{:}\n{:}\n{:}\n{:} -99999 -99999 -99999 -99999\n1\n3\n"
                    .format(pick_onset, pick.phase_hint, pick_polarity,
                            pick_weight))
            w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> "
                    "<NOT_SET> {:} <NOT_SET>\n".format(pick_create_time))
            w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> {:} <NOT_SET>\n"
                    .format(pick_create_time))
        w.write("[-PICKSET]\n")

        # Insert event
        w.write("[+EVENTSET]\n")
        w.write("event\n")
        orig_time = (str(origin.time.datetime).replace(":", "")
                     .replace("-", "").replace(".", "").replace(" ", "")[:-3])
        orig_time_microsec = (str(pick.time.datetime).replace(":", "")
                              .replace("-", "").replace(" ", "_")[:-3])
        w.write("{:}\n".format(orig_time))
        w.write("{:}_{:}_{:}\n".format(
            orig_time_microsec, int(origin.latitude), int(origin.longitude)))
        w.write("<NOT_SET>\n<NOT_SET>\n18\n2\n<NOT_SET> 7\n")
        origin_createtime = (str(origin.creation_info.creation_time)
                             .replace("-", ".").replace(" ", "-")
                             .replace("T",                   "-")[:-3])
        w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> "
                "{:} <NOT_SET>\n".format(origin_createtime))
        w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> {:} <NOT_SET>\n"
                .format(origin_createtime))
        w.write("#00661d\n{:}\n".format(len(picked_stations)))
        for pick_id in picked_stations:
            w.write("{:}\n".format(pick_id))
        w.write("[-EVENTSET]\n")
        w.write("[+ORIGINSET]\norigin\n")
        w.write("{:}_{:}_{:}\n".format(
            orig_time_microsec, int(origin.latitude), int(origin.longitude)))
        origin_time = (str(event.origins[0].time)
                       .replace("-", ".").replace("T", "-")[:-1]
                       .replace(" ", "-")[:-3])
        w.write("{:} -99999 -99999 -99999 -99999\n".format(origin_time))
        w.write("{:} -99999 -99999 -99999 -99999\n".format(origin.latitude))
        w.write("{:} -99999 -99999 -99999 -99999\n".format(origin.longitude))
        w.write("{:} -99999 -99999 -99999 -99999\n".format(origin.depth/1000))
        w.write("5\nF\nF\n<NOT_SET>\nhypo71\nlocal\n")
        w.write("-99999 -99999 -99999 -99999 -99999 -99999 -99999 -99999 "
                "-99999 -99999 -99999 -99999 -99999 -99999 -99999 -99999 "
                "-99999 -99999 -99999 -99999 -99999 -99999 -99999 -99999 "
                "-99999 -99999 -99999\n")
        Np = [arrival for arrival in origin.arrivals if arrival.phase == "P"]
        w.write("{0:} {0:} {1:} {1:} -99999 0 {2:} -99999 "
                "<NOT_SET> {3:} {4:} {5:}\n".format(
                    origin.quality.used_phase_count, len(Np),
                    origin.quality.azimuthal_gap,
                    origin.quality.minimum_distance,
                    origin.quality.maximum_distance,
                    origin.quality.median_distance))
        w.write("6\n0\n6\n")
        w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> "
                "{:} <NOT_SET>\n".format(origin_createtime))
        w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> {:} <NOT_SET>\n"
                .format(origin_createtime))
        w.write("[-ORIGINSET]\n")

        # Insert stations
        w.write("[+STATIONSET]\n")
        stations = open(station_file, "r")
        for station in stations:
            w.write("station\n")
            w.write("{:}_20000101_000000_20241231_000000\n"
                    .format(station.split()[0]))
            w.write("{:}\n".format(station.split()[0]))
            w.write("{:7.4f} -99999 -99999 -99999 -99999\n"
                    .format(float(station.split()[1])))
            w.write("{:7.4f} -99999 -99999 -99999 -99999\n"
                    .format(float(station.split()[2])))
            w.write("{:4g} -99999 -99999 -99999 -99999\n"
                    .format(float(station.split()[3])))
            w.write("2007.01.01-00:00:00.000 -99999 -99999 -99999 -99999\n")
            w.write("2024.12.31-00:00:00.000 -99999 -99999 -99999 -99999\n")
            w.write("<NOT_SET>\n")
            w.write("<NOT_SET>\n")
            w.write("VOILA%0A <NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> "
                    "<NOT_SET> 2018.01.11-13:28:30.021 <NOT_SET>\n")
            w.write("<NOT_SET> <NOT_SET> <NOT_SET> <NOT_SET> "
                    "2018.01.11-13:28:30.021 <NOT_SET>\n")
        w.write("[-STATIONSET]\n")
        stations.close()
        w.write("end")
