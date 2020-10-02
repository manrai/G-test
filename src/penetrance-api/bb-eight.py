# Arjun Manrai
# NIH/NHLBI BDC
# BB-EIGHT Penetrance API

from flask import Flask, request, jsonify, Response
import numpy as np
import pandas as pd
import sqlite3
import json
import subprocess

app = Flask(__name__)

# Establish SQL connection to DB
conn_1 = sqlite3.connect('../../data/prebuilt-dbs/variants.db')
conn_2 = sqlite3.connect('../../data/prebuilt-dbs/bdc-test.db')
df_1 = pd.read_sql_query("SELECT * from variants", conn_1) # all variants
df_2 = pd.read_sql_query("SELECT * from clinvar", conn_1) # all bdc
df_3 = pd.read_sql_query("SELECT * from bdc", conn_2) # all bdc

# No parameters specified, return info in databases and help message
@app.route('/', methods=['GET'])
def get_welcome():
    msg = '''
            <img src="https://biodatacatalyst.nhlbi.nih.gov/static/bdc-logo-45457c10c57bf2d4eba026b0d79b61f9.svg">
            <h1>Welcome to the NHLBI BioData Catalyst BB-EIGHT API</h1>
            <p>This API allows investigators to specify phenotype and genotype
            criteria dynamically to obtain penetrance estimates for genetic variants
            across populations.</p>
            For more information please visit:
            <a href="https://github.com/manrai/bb-eight">BB-EIGHT codebase</a>
            '''
    return msg

# Return variant by rsID
@app.route('/variants/rsid/<string:rsid>', methods=['GET'])
def get_variant_by_rsid(rsid):
    variants = []
    for index, row in df_1.iterrows():
        if row["rsID"] == rsid:
            r = {
                'Chromosome': row["Chromosome"],
                'Position': row["Position"],
                'rsID': row["rsID"],
                'Reference': row["Reference"],
                'Alternate': row["Alternate"],
                'Frequency': row["Allele.Frequency"]
            }
            variants += [r]
    return jsonify(variants)

# Return variant by chr,pos,ref,alt
@app.route('/variants/position', methods=['GET'])
def get_variant_by_position():
    chr = int(request.args.get('chr',None))
    pos = int(request.args.get('pos',None))
    ref = str(request.args.get('ref',None))
    alt = str(request.args.get('alt',None))
    variants = []
    for index, row in df_1.iterrows():
        if row["Chromosome"] == chr and row["Position"] == pos and \
        row["Reference"] == ref and row["Alternate"] == alt:
            r = {
                'Chromosome': row["Chromosome"],
                'Position': row["Position"],
                'rsID': row["rsID"],
                'Reference': row["Reference"],
                'Alternate': row["Alternate"],
                'Frequency': row["Allele.Frequency"]
            }
            variants += [r]
    return jsonify(variants)

# Return random variant
@app.route('/variants/random', methods=['GET'])
def get_random_variant():
    row = df_1.sample()
    variants = []
    r = {
        'Chromosome': int(row["Chromosome"].values[0]),
        'Position': int(row["Position"].values[0]),
        'rsID': row["rsID"].values[0],
        'Reference': row["Reference"].values[0],
        'Alternate': row["Alternate"].values[0],
        'Frequency': float(row["Allele.Frequency"].values[0])
    }
    return jsonify(r)

# Return all variants
@app.route('/variants/all', methods=['GET'])
def get_all_variants():
    variants = []
    for index, row in df_1.iterrows():
        r = {
            'Chromosome': row["Chromosome"],
            'Position': row["Position"],
            'rsID': row["rsID"],
            'Reference': row["Reference"],
            'Alternate': row["Alternate"],
            'Frequency': row["Allele.Frequency"]
        }
        variants += [r]
    return jsonify(variants)

# Return clinvar matches by position
@app.route('/clinvar/position', methods=['GET'])
def get_clinvar_by_position():
    chr = int(request.args.get('chr',None))
    pos = int(request.args.get('pos',None))
    variants = []
    for index, row in df_2.iterrows(): # df_2 = clinvar
        if row["Chromosome"] == chr and row["Position"] == pos:
            r = {
                'Chromosome': row["Chromosome"],
                'Position': row["Position"],
                'Reference': row["Reference"],
                'Alternate': row["Alternate"],
                'Pathogenic': row["Pathogenic"]
            }
            variants += [r]
    return jsonify(variants)

# return penetrance MAP, low, high estimate from user G & P criteria
@app.route('/penetrance', methods=['POST'])
def compute_penetrance():
    # parameters
    # vus, likely_pathogenic, pathogenic [g]
    # lvwt_min, lvwt_max, age_min, age_max [p]
    # gender, race [p]
    # cohorts [g, p]
    # htn [p]

    # get parameters
    req_data = request.get_json()

    # define which variant classes are considered in pen. calculation
    vus = 0 # default
    likely_pathogenic, pathogenic = 1, 1 # default
    if 'vus' in req_data: vus = req_data['vus']
    if 'likely_pathogenic' in req_data: likely_pathogenic = req_data['likely_pathogenic']
    if 'pathogenic' in req_data: pathogenic = req_data['pathogenic']

    # phenotype criteria for CONTROLS
    lvwt_min, lvwt_max = 8, 13 # default
    age_min, age_max = 18, 40 # default
    if 'lvwt_min' in req_data: lvwt_min = req_data['lvwt_min']
    if 'lvwt_max' in req_data: lvwt_max = req_data['lvwt_max']
    if 'age_min' in req_data: age_min = req_data['age_min']
    if 'age_max' in req_data: age_max = req_data['age_max']

    # demo. for CONTROLS
    gender, race = "all", "all"
    if 'gender' in req_data: gender = req_data['gender']
    if 'race' in req_data: race = req_data['race']

    # cohorts for CONTROLS
    cohorts = ["jackson","cardia"]
    if 'cohorts' in req_data: cohorts = req_data['cohorts']

    # htn for CONTROLS
    htn = 1
    if 'htn' in req_data: htn = req_data['htn']

    # call compute_penetrance.R
    command = "Rscript compute_penetrance.R"
    command += " --vus " + str(vus)
    command += " --likely_pathogenic " + str(likely_pathogenic)
    command += " --pathogenic " + str(pathogenic)

    out = subprocess.run(command,shell=True,capture_output=True)
    vals = str(out.stdout)
    vals = vals.split('"')[1].split(';')

    low = float(vals[0])
    map = float(vals[1])
    high = float(vals[2])

    s_1 = {
        "Low": low,
        "MAP": map,
        "High": high
    }
    response = Response(json.dumps(s_1), status = 200, mimetype="application/json")
    return response

if __name__ == '__main__':
    app.run(port=5000)
