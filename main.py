"""Generate the figures and static webpages for the long term LCM validation

Copyright William Clarke, University of Oxford, 2022
"""

from pathlib import Path
import json
from distutils.version import LooseVersion

import pandas as pd
from jinja2 import Environment, FileSystemLoader

from single_version_plots import generate_version_plots
import timecourse_plots

# 1. ID all the folders that relate to packages/algorithms/models
res_dir = Path('results')
pkg_list = []
for fp in res_dir.rglob('pkg_info.json'):
    with open(fp) as pkinfo:
        current_dict = json.load(pkinfo)
    current_dict['location'] = str(fp.parent)
    current_dict['name'] = str(fp.parent.stem)
    ver_list = [x.stem.replace('_', '.') for x in fp.parent.glob('*.gz')]
    ver_list.sort(key=LooseVersion)
    current_dict['versions'] = ver_list
    current_dict['current_version'] = ver_list[-1]
    pkg_list.append(current_dict)

# 2. Generate figures for each
# 2a. Any new version figures
html_path = Path('html')
out_dir = html_path / 'figures'
out_dir.mkdir(exist_ok=True, parents=True)

for pkg in pkg_list:
    latest_version = pkg['versions'][-1].replace('.', '_')
    file_loc = Path(pkg['location']) / (latest_version + '.gz')
    answer_df = pd.read_csv('answers.gz', index_col=0)
    answer_df.columns = [str(i) for i in range(1, 22)]
    res_df = pd.read_csv(file_loc, index_col=0)
    name = pkg['name']
    file_outputs = generate_version_plots(res_df, answer_df, name, out_dir)
    pkg.update(file_outputs)

# 2b. New timeline figure
answer_df = pd.read_csv('answers.gz', index_col=0).T
answer_df.index = [str(i) for i in range(1, 22)]
answer_df.drop(['lipids', 'Ace'], axis=1, inplace=True)
# print(res_df)
for pkg in pkg_list:
    file_loc = Path(pkg['location'])
    curr_metrics = []
    for ver in pkg['versions']:
        res_df = pd.read_csv(
            file_loc / (ver.replace('.', '_') + '.gz'),
            index_col=0).T
        res_df.drop(['Mac'], axis=1, inplace=True)
        curr_metrics.append(
            timecourse_plots.calculate_metrics(
                res_df,
                answer_df))

    name = pkg['name']
    file_outputs = timecourse_plots.generate_tc_plots(
        curr_metrics,
        pkg['versions'],
        name,
        out_dir)
    pkg.update(file_outputs)


# 3. Generate HTML output
environment = Environment(loader=FileSystemLoader("templates/"))

# 3a. ID any algos that have new versions or are new and (re)generate any HTML
template = environment.get_template("pkg_results.html")
for pkg in pkg_list:
    content = template.render(pkg)
    html_link = f'{pkg["name"]}.html'
    pkg['link'] = str(html_link)
    with open(html_path / html_link, mode="w", encoding="utf-8") as message:
        message.write(content)

# 3b. Generate landing page.
template = environment.get_template("index.html")
content = template.render(packages=pkg_list)
with open(html_path / 'index.html', mode="w", encoding="utf-8") as message:
    message.write(content)
