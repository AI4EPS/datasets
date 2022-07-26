# %%
import urllib
import urllib.request as request
import re, os
from glob import glob
from tqdm import tqdm
import collections

# %%
root_dir = "catalogs"
root_url = "http://ncedc.org/ftp/pub/catalogs/ncss/hypoinverse/phase2k"
if not os.path.exists(root_dir):
    os.mkdir(root_dir)

# %%
def get_years(root_url):
    html = urllib.request.urlopen(root_url).read().decode()
    pattern = re.compile("<a href=\"\d\d\d\d/\">", re.S)
    tmp_years = re.findall(pattern, html)
    years = [re.findall("\d\d\d\d", yr)[0] for yr in tmp_years]
    year_urls = {yr: root_url+"/"+yr for yr in years}
    return year_urls

year_urls = get_years(root_url)

# %%
def get_files(year_urls):
    file_urls = {}
    for year, url in year_urls.items():
        html = urllib.request.urlopen(url).read().decode()
        pattern = re.compile("<a href=\".*?\.phase\.Z\">", re.S)
        tmp_files = re.findall(pattern, html)
        files = [re.findall("\d.*?\.phase\.Z", fl)[0] for fl in tmp_files]
        file_urls[year] = [url+"/"+fl for fl in files]
    return file_urls

file_urls = get_files(year_urls)
# %%
def download_files(file_urls, root_dir):
    for year in file_urls:
        data_dir = os.path.join(root_dir, year)
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        for url in file_urls[year]:
            print("Downloading: "+url)
            request.urlretrieve(url, os.path.join(data_dir, url.split('/')[-1]))
            os.system("uncompress "+os.path.join(data_dir, url.split('/')[-1]))

download_files(file_urls, root_dir)
        
# %%
def merge_files(file_urls, root_dir, fout):

    catlog = []
    for year in file_urls:
        for url in file_urls[year]:
            print("Merging {}".format(url))
            with open(os.path.join(root_dir, year, url.split('/')[-1].rstrip(".Z")), 'r') as fp:
                lines = fp.readlines()
                catlog += lines

    with open(fout, 'w') as fp:
      for line in catlog:
        fp.write(line)
    print(f"Finish writing {len(catlog)} lines to {fout}")

merge_files(file_urls, root_dir, "catalogs.txt")

# %%
def build_dict(catalog):
    dataset1 = collections.OrderedDict()
    with open(catalog) as fp:
        for line in fp:
            if line[0].isspace():
                continue
            elif len(line) > 130:
                event_id = line
                dataset1[event_id] = []
            elif len(line) <= 130:
                dataset1[event_id].append(line)
            else:
                print("Unrecognized line: %s" % line)

    # dataset organized by event then station
    dataset2 = collections.OrderedDict()
    for event in tqdm(dataset1, desc="build_dict: "):
        stationset = collections.OrderedDict()
        for line in dataset1[event]:
            # if line[111:113] not in ["  ", "--", "00"]:
            #     sta_id = line[:7]+line[111:113]#plus location code
            # else:
            #     sta_id = line[:7]+"--"
            sta_id = line[:7]+line[111:113]
            if sta_id in stationset:
                stationset[sta_id].append(line)
            else:
                stationset[sta_id] = [line]
        dataset2[event] = stationset

    return dataset2

catalog_dict = build_dict("catalogs.txt")
# %%
def extract_ps(catalog):

    # pick the best P and S pickers
    dataset = collections.OrderedDict()
    for event in tqdm(catalog, desc="extract P/S picks:"):
        stationset = collections.OrderedDict()
        for sta_id in catalog[event]:
            best_p = 10
            best_s = 10
            id_p = 0
            id_s = 0
            found_p = 0
            found_s = 0
            for j, line in enumerate(catalog[event][sta_id]):
                if line[14] == 'P':
                    if int(line[16]) < best_p:
                        best_p = int(line[16])
                        id_p = j
                        found_p = 1
                if line[47] == 'S':
                    if int(line[49]) < best_s:
                        best_s = int(line[20])
                        id_s = j
                        found_s = 1

            if found_p and found_s:
                stationset[sta_id] = [catalog[event][sta_id][id_p], catalog[event][sta_id][id_s]]

            dataset[event] = stationset
    
    return dataset

dataset_ps = extract_ps(catalog_dict)

# %%
def write_ps(dataset, fname, with_event=False):
    with open(fname, 'w') as fp:
        for event in tqdm(dataset, desc="write_data:"):
            if with_event and len(dataset[event]) > 0:
                fp.write(event)
            for sta_id in dataset[event]:
                for line in dataset[event][sta_id]:
                    fp.write(line)

write_ps(dataset_ps, "catalogs_ps.txt", with_event=True)
# %%
