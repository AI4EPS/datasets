# %%
from re import S
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import obspy
from collections import namedtuple, OrderedDict
from obspy.clients.fdsn import Client
client = Client("NCEDC")
import pandas as pd

#%%
def read_ps_catalog(fname):
    events = OrderedDict()
    phases = OrderedDict()
    index = -1
    with open(fname) as fp:
        for line in tqdm(fp, desc="read catalog"):
            if len(line) > 130:
                index += 1
                event_line = line
                events[index] = event_line
                phases[index] = []
            elif len(line) <= 130:
                phase_line = line
                phases[index].append(phase_line)
            else:
                print("Unrecognized line: %s" % line)
    return events, phases

file_catalog = "catalogs_ps.txt"
events, phases = read_ps_catalog(file_catalog)

# %%
def to_float(string):
  if string.strip() == '':
    return 0
  else:
    return float(string.strip())

Event = namedtuple("event", ["time", "latitude", "longitude",
                             "depth_km", "magnitude", "magnitude_type", "index"])
def read_event_line(line, idx):
    time = line[:4]+"-"+line[4:6]+"-"+line[6:8]+"T"+line[8:10]+":"+line[10:12]+":"+line[12:14]+"."+line[14:16] 
    latitude = to_float(line[16:18]) + to_float(line[19:23])/6000.0
    longitude = -(to_float(line[23:26]) + to_float(line[27:31])/6000.0)
    if line[18] == 'S':
      latitude *= -1.0
    if line[26] == 'E':
      longitude *= -1.0
    depth_km = to_float(line[31:36])/100.0
    magnitude_type = (line[146:147]).strip()
    magnitude = to_float(line[147:150])/100.0
    return Event(time=time, latitude=np.round(latitude,4), longitude=np.round(longitude,4), 
                 depth_km=depth_km, magnitude=magnitude, magnitude_type=magnitude_type, index=idx)

def read_p_pick(line):
    if float(line[30:32].replace(" ", "0")) < 60:
      tp = (line[17:21]+"-"+line[21:23]+"-"
            +line[23:25]+"T"+line[25:27]+":"
            +line[27:29]+":"+line[30:32]+'.'
            +line[32:34]).replace(" ", "0")
      tp = obspy.UTCDateTime(tp)
    else:
      tp = (line[17:21]+"-"+line[21:23]+"-"
            +line[23:25]+"T"+line[25:27]+":"
            +line[27:29]+":"+'00'+'.'
            +line[32:34]).replace(" ", "0")
      tp = obspy.UTCDateTime(tp) + float(line[30:32].replace(" ", "0"))
    remark = line[13:15].strip()
    weight = line[16:17].strip()
    channel = line[9:12].strip()
    return tp, remark, weight, channel

def read_s_pick(line):
    if float(line[42:44].replace(" ", "0")) < 60:
      ts = (line[17:21]+"-"+line[21:23]+"-"
            +line[23:25]+"T"+line[25:27]+":"
            +line[27:29]+":"+line[42:44]+'.'
            +line[44:46]).replace(" ", "0")
      ts = obspy.UTCDateTime(ts)
    else:
      ts = (line[17:21]+"-"+line[21:23]+"-"
            +line[23:25]+"T"+line[25:27]+":"
            +line[27:29]+":"+"00"+'.'
            +line[44:46]).replace(" ", "0")
      ts = obspy.UTCDateTime(ts) + float(line[42:44].replace(" ", "0"))
    remark = line[46:48].strip()
    weight = line[49:50].strip()
    channel = line[9:12].strip()
    return ts, remark, weight, channel

Pick = namedtuple("pick",["p_time", "p_remark", "p_weight", "p_channel",
                          "s_time", "s_remark", "s_weight", "s_channel",
                          "first_motion", "distance_km",
                          "emergence_angle", "azimuth", 
                          "network", "station", "location_code", "event_index"])
def read_phase_line(p_line, s_line, index):

    line = p_line
    network = (line[5:7]).strip()
    station = (line[:5]).strip()
    location_code = (line[111:113]).strip()
    distance_km = to_float(line[74:78])/10.0
    emergence_angle = to_float(line[78:81])
    # duration = to_float(line[87:91])
    azimuth = to_float(line[91:94])
    first_motion = (line[15:16]).strip()
    p_time, p_remark, p_weight, p_channel = read_p_pick(p_line)
    s_time, s_remark, s_weight, s_channel = read_s_pick(s_line)
    
    return Pick(p_time=p_time, p_remark=p_remark, p_weight=p_weight, p_channel=p_channel,
                s_time=s_time, s_remark=s_remark, s_weight=s_weight, s_channel=s_channel,
                first_motion=first_motion, distance_km=distance_km, 
                emergence_angle=emergence_angle, azimuth=azimuth, 
                network=network, station=station, location_code=location_code, event_index=index)

def resample(stream, default_sampling_rate):
    sampling_rate = stream[-1].stats.sampling_rate
    if sampling_rate < default_sampling_rate:
        print("Resample %s" % stream)
        stream.resample(default_sampling_rate) ## resample to 100HZ
        print("After resample: %s" % stream)
    elif sampling_rate > default_sampling_rate:
        print("Resample %s" % stream)
        if np.mod(sampling_rate, default_sampling_rate) == 0: ##directly throw away the data
            stream.decimate(int(sampling_rate//default_sampling_rate), strict_length=False, no_filter=True)
        else:
            stream.resample(default_sampling_rate) ## resample to 100HZ
        print("After resample: %s" % stream)
    return stream

Station = namedtuple("station", ["id", "latitude", "longitude", "elevation_m", "unit"])
def download_waveform(pick, waveform_path):

    Tstart =  pick.p_time - 60.0
    Tend = pick.p_time + 60.0
    
    if pick.p_channel[:-1] != pick.s_channel[:-1]:
        channels = [pick.p_channel[:-1], pick.s_channel[:-1]]
    else: 
        channels = [pick.p_channel[:-1]]
    channels = set(channels + ["HN", "HH", "EH", "BH", "DP"])

    stream_list = []
    station_list = []
    for channel in channels:
        fname = pick.network+"."+pick.station+"."+pick.location_code+"."+channel+"."+f"{pick.event_index:07d}"+".mseed"
        # if os.path.exists(os.path.join(waveform_path, fname)):
        #     stream = obspy.read(os.path.join(waveform_path, fname))
        #     stream_list.append(stream)
        # else:
        try:
            stream = client.get_waveforms(network=pick.network, station=pick.station, location=pick.location_code, 
                                        channel=channel+"?", starttime=Tstart, endtime=Tend)
            stream = stream.detrend('linear')
            stream = stream.merge(fill_value=0)
            stream = stream.trim(Tstart, Tend, pad=True, fill_value=0)
            stream = stream.sort()

            station = client.get_stations(network=pick.network, station=pick.station, location=pick.location_code, 
                                        channel=channel+"?", starttime=Tstart, endtime=Tend, level="response")
            stream.attach_response(station)
            stream.remove_sensitivity()
            coord = station.get_coordinates(stream[0].get_id(), datetime=Tstart)
            response = station.get_response(stream[0].get_id(), datetime=Tstart)
            sta_id = pick.network+"."+pick.station+"."+pick.location_code+"."+channel
            unit = response.instrument_sensitivity.input_units.lower()
            station = Station(id=sta_id,
                                latitude=np.round(coord["latitude"],4), longitude=np.round(coord["longitude"],4), 
                                elevation_m=coord["elevation"], unit=unit)
            
            stream = resample(stream, 100) ## resample to 100 Hz
            station_list.append(station)
            stream_list.append(stream)
            stream.write(os.path.join(waveform_path, fname),  format="MSEED")

        except Exception as e:
            if str(e)[:len("No data available")] != "No data available":
                print("Failed downloading: "+fname, "Error: "+str(e))

    return stream_list, station_list

def calc_snr(vec, anchor, dt):
    npts = int(3/dt)
    eps = 10*np.finfo(vec.dtype).eps
    snr = (np.std(vec[anchor:anchor+npts, :], axis=0)+eps) / (np.std(vec[anchor-npts:anchor, :], axis=0)+eps)
    return snr

def data2vec(i, data, vec, shift, window_size):
    if shift >= 0:
        if len(data)+shift <= window_size:
            vec[shift:len(data)+shift, i] = data
        else:
            vec[shift:, i] = data[:window_size-shift]
    else:
        if len(data)+shift <= window_size:
            vec[:len(data)+shift, i] = data[-shift:]
        else:
            vec[:, i] = data[-shift:window_size-shift]
    return vec

Extra = namedtuple("pick_extra", ["p_idx", "s_idx", "channels", "snr", "dt", "station_index", "latitude", "longitude", "elevation_m", "unit"])
def convert_sample(pick, event, stream, station, sample_path):

    dt = stream[-1].stats.delta
    npts = stream[-1].stats.npts 
    starttime = stream[-1].stats.starttime                                  
    endtime = stream[-1].stats.endtime

    p_idx = int(np.around( (pick.p_time - starttime)/(endtime - starttime)*npts )) 
    s_idx = int(np.around( (pick.s_time - starttime)/(endtime - starttime)*npts ))
    
    anchor = 6000
    window_size = 12000
    vec = np.zeros([window_size, 3])
    shift = anchor - p_idx
    if np.abs(shift) <= 1:
        shift = 0
    p_idx += shift
    s_idx += shift

    order = ['3','2','1','E','N','Z']
    order = {key: i for i, key in enumerate(order)}
    comps = [x.get_id() for x in stream]
    comps = sorted(comps, key=lambda x: order[x[-1]])

    fname = comps[0][:-1]+"."+f"{pick.event_index:07d}"+".npz"
    if len(comps) == 3:
        for i, c in enumerate(comps):
            data = stream.select(id=c)[0].data
            data2vec(i, data, vec, shift, window_size)
    elif len(comps) < 3:
        for c in comps:
            if c[-1] == "E":
                i = 0
            elif c[-1] == "N":
                i = 1
            elif c[-1] == "Z":
                i = 2
            else:
                print(f"Unknown channels: {comps}")
                return (1, None)
            data = stream.select(id=c)[0].data
            data2vec(i, data, vec, shift, window_size)
    else:
        print(f"Unknown channels: {comps}")
        return (1, None)

    snr = calc_snr(vec, anchor, dt)
    channels = ",".join([x.split(".")[-1] for x in comps])
    extra = Extra(p_idx=p_idx, s_idx=s_idx, snr=tuple(snr.tolist()), dt=dt,
                  station_index=station.id, channels=channels,
                  latitude=station.latitude, longitude=station.longitude, elevation_m=station.elevation_m,
                  unit=station.unit)
    np.savez(os.path.join(sample_path, fname), 
            data=vec.astype("float32"), dt=dt, p_idx=p_idx, s_idx=s_idx, snr=snr.tolist(), 
            p_time=pick.p_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3], s_time=pick.s_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3],
            p_remark=pick.p_remark, p_weight=pick.p_weight, s_remark=pick.s_remark, s_weight=pick.s_weight,
            first_motion=pick.first_motion, distance_km=pick.distance_km, 
            azimuth=pick.azimuth, emergence_angle=pick.emergence_angle, 
            network=pick.network, station=pick.station, location_code=pick.location_code, 
            station_latitude=station.latitude, station_longitude=station.longitude, station_elevation_m=station.elevation_m,
            event_latitude=event.latitude, event_longitude=event.longitude, event_elevation_m=event.depth_km,
            event_time=event.time, event_magnitude=event.magnitude, event_magnitude_type=event.magnitude_type,
            unit=station.unit, channels=channels, event_index=pick.event_index)
            
    return (0, extra)

# %%
def create_dataset(idx, events, phases, waveform_path, sample_path):
    events_ = []
    phases_ = []
    extras_ = []
    for i in tqdm(idx, desc="create dataset"):
        event = read_event_line(events[i], i)
        # events_.append(event)
        phase_lines = phases[i]
        has_phase = False
        for j in range(0, len(phase_lines), 2):
            pick = read_phase_line(phase_lines[j], phase_lines[j+1], i)
            # phases_.append(pick)
            (waveforms, stations) = download_waveform(pick, waveform_path)
            for w, s in zip(waveforms, stations):
                (status, extra) = convert_sample(pick, event, w, s, sample_path)
                if (status == 0):
                    phases_.append(pick)
                    extras_.append(extra)
                    has_phase = True
        if has_phase:
            events_.append(event)

    return events_, phases_, extras_

waveform_path = "waveforms"
sample_path = "dataset"
if not os.path.exists(waveform_path):
    os.mkdir(waveform_path)
if not os.path.exists(sample_path):
    os.mkdir(sample_path)
events_tuple, phases_tuple, extra_tuple = create_dataset([len(events)-1], events, phases, waveform_path, sample_path)
# events_tuple, phases_tuple, extra_tuple = create_dataset(range(len(events)//6, len(events)), events, phases, waveform_path, sample_path)

if len(events_tuple) >= 1:
    events_df = pd.DataFrame(data=events_tuple)
    events_df = events_df.set_index("index")
    events_df.to_csv("event_catalog.csv", sep="\t")

if len(phases_tuple) >= 1:
    phases_df = pd.DataFrame(data=phases_tuple)
    extra_df = pd.DataFrame(data=extra_tuple)
    phases_df["p_idx"] = extra_df["p_idx"]
    phases_df["s_idx"] = extra_df["s_idx"]
    phases_df["dt"] = extra_df["dt"]
    phases_df["latitude"] = extra_df["latitude"]
    phases_df["longitude"] = extra_df["longitude"]
    phases_df["elevation_m"] = extra_df["elevation_m"]
    phases_df["channels"] = extra_df["channels"]
    phases_df["unit"] = extra_df["unit"]
    phases_df["snr"] = extra_df["snr"].apply(lambda x: ",".join([f"{i:.2f}" for i in x]))
    phases_df["p_time"] = phases_df["p_time"].apply(lambda x: x.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3])
    phases_df["s_time"] = phases_df["s_time"].apply(lambda x: x.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3])
    phases_df.to_csv("phases.csv", sep="\t", index=False,
        columns=["network","station","location_code", 
        "p_idx", "p_time","p_remark","p_weight",
        "s_idx", "s_time","s_remark","s_weight",
        "first_motion","distance_km","emergence_angle","azimuth",
        "latitude", "longitude", "elevation_m", "unit", "dt",
        "event_index", "channels", "snr"])

# %%
import glob
def check_sample(sample_path):
    files = glob.glob(sample_path+"/*.npz")
    for f in files:
        meta = np.load(f)
        for k in meta.files:
            print(k, meta[k])
        print(meta["data"].dtype)
        
        plt.figure()
        plt.plot(meta["data"][:,-1])
        ylim = plt.ylim()
        plt.plot([meta["p_idx"], meta["p_idx"]], ylim)
        plt.plot([meta["s_idx"], meta["s_idx"]], ylim)
        plt.show()

        return

sample_path = "dataset"
check_sample(sample_path)

# %%

# %%
