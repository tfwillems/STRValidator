#!/usr/bin/env python2.7

from flask import Flask, url_for, make_response, render_template, request
from flask import json, jsonify, current_app, Blueprint
from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy.sql import func

from json import dumps
from sqlalchemy.orm import class_mapper
from sqlalchemy import or_, and_, ForeignKey
from sqlalchemy import func

from flask.ext.compress import Compress
from flask.ext.assets import Environment, Bundle

from subprocess import Popen,PIPE

import argparse
import os.path
import collections

admin = Blueprint('admin', __name__, static_folder='static', url_prefix='/')
db    = SQLAlchemy()

class HaploidCalls(db.Model):
    __tablename__ = 'haploid_calls'
    id         = db.Column(db.Integer, primary_key=True)
    chrom      = db.Column(db.Text)
    start      = db.Column(db.Integer)
    end        = db.Column(db.Integer)
    total_one  = db.Column(db.Integer)
    total_two  = db.Column(db.Integer)
    g1_1       = db.Column(db.Integer)
    g1_2       = db.Column(db.Integer)
    g2_1       = db.Column(db.Integer)
    g2_2       = db.Column(db.Integer)
    sample_1   = db.Column(db.Text)
    sample_2   = db.Column(db.Text)

class DiploidCalls(db.Model):
    __tablename__ = 'diploid_calls'
    id         = db.Column(db.Integer, primary_key=True)
    chrom      = db.Column(db.Text)
    start      = db.Column(db.Integer)
    end        = db.Column(db.Integer)
    total_one  = db.Column(db.Integer)
    total_two  = db.Column(db.Integer)
    g1_1       = db.Column(db.Integer)
    g1_2       = db.Column(db.Integer)
    g2_1       = db.Column(db.Integer)
    g2_2       = db.Column(db.Integer)
    sample_1   = db.Column(db.Text)
    sample_2   = db.Column(db.Text)

@admin.route('/')
def index():
    return render_template("index.html", comparisons=current_app.config['COMPARISONS'].keys(), analyses=current_app.config['ANALYSES'])

@admin.route('getbubbleplot')
def get_bubbleplot():
    return render_template("bubbleplot.html")

@admin.route('getbarplot')
def get_barplot():
    return render_template("barplot.html", min_bp_diff=current_app.config['MIN_BP_DIFF'], max_bp_diff=current_app.config['MAX_BP_DIFF'])
    
@admin.route('getplot/<plot_type>')
def get_plot(plot_type):
    if plot_type == "Bubble plot":
        return get_bubbleplot()
    elif plot_type == "Bar plot":
        return get_barplot()
    else:
        exit("Invalid plot type %s"%(plot_type))

def get_alignment_info(chroms, starts, stops, sample_lists, annot_lists):
    input        = "\n".join(map(lambda x: "\t".join(x), zip(chroms, starts, stops, sample_lists, annot_lists)))
    cmd          = "%s --bams %s --indexes %s --out /dev/stdout --regions /dev/stdin --fasta %s"%(current_app.config['vizalign'], current_app.config['bams'],
                                                                                                  current_app.config['bais'], current_app.config['fasta_dir'])
    stream       = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE)
    res          = stream.communicate(input=input)[0]
    return res

@admin.route('getbubbleinfo')
def get_bubble_info():
    schema = current_app.config['ACTIVE_DB']
    res    = db.session.query(schema.total_one, schema.total_two, func.count(schema.total_one)).group_by(schema.total_one, schema.total_two).all()
    res    = sorted(res, key=lambda x:-x[2])
    return jsonify(result=res)

@admin.route('getdiffinfo')
def get_diff_info():
    schema      = current_app.config['ACTIVE_DB']
    res         = db.session.query(schema.total_one, schema.total_two, func.count(schema.total_one)).group_by(schema.total_one, schema.total_two).all()
    diff_counts = collections.defaultdict(int)
    for item in res:
        diff_counts[item[0]-item[1]] += item[2]
    return jsonify(result=diff_counts.items())

def extract_alignments(query):
    chroms, starts, stops, sample_lists, annot_lists = [], [], [], [], []
    for record in query:
        chroms.append(record.chrom)
        starts.append(str(record.start))
        stops.append(str(record.end))
        if record.sample_1 == record.sample_2:
            sample_lists.append(record.sample_1)
            annot_lists.append(str(record.g1_1)+"/"+str(record.g1_2))
        else:
            sample_lists.append(record.sample_1 + "," + record.sample_2)
            annot_lists.append(str(record.g1_1)+"/"+str(record.g1_2)+ "," + str(record.g2_1)+"/"+str(record.g2_2))

    res      = get_alignment_info(chroms, starts, stops, sample_lists, annot_lists)
    response = make_response(res)
    response.headers.add('Access-Control-Allow-Origin','*')
    return response


@admin.route('getbubblealignments/<x>/<y>')
def get_bubble_alignments(x, y):
    x, y   = int(x), int(y)
    schema = current_app.config['ACTIVE_DB']
    res = db.session.query(schema).filter(schema.total_one == x, schema.total_two == y)
    return extract_alignments(res)

@admin.route('getdiffalignments/<diff>')
def get_diff_alignments(diff):
    diff   = int(diff)
    schema = current_app.config['ACTIVE_DB']
    if diff == current_app.config['MIN_BP_DIFF']:
        res = db.session.query(schema).filter(schema.total_one - schema.total_two <= diff)
    elif diff == current_app.config['MAX_BP_DIFF']:
        res = db.session.query(schema).filter(schema.total_one - schema.total_two >= diff)
    else:
        res = db.session.query(schema).filter(schema.total_one - schema.total_two == diff)
    return extract_alignments(res)

@admin.route('set_comparison/<comp>')
def set_comparison(comp):
    current_app.config['ACTIVE_DB'] = current_app.config['COMPARISONS'][comp]
    return jsonify(success="TRUE")


def create_app(bams, bais, fasta_dir, vizalign):
    app = Flask(__name__)

    # Check that .bam files exist
    bams = bams.split(",")
    for file in bams:
        if not os.path.isfile(file):
            exit("ERROR: BAM file %s does not exist"%(file))
    app.config['bams'] = ",".join(bams)

    # Check that .bai files exist
    bais = bais.split(",")
    for file in bais:
        if not os.path.isfile(file):
            exit("ERROR: BAI file %s does not exist"%(file))
    app.config['bais'] = ",".join(bais)

    # Check that fasta directory exists
    if not os.path.isdir(fasta_dir):
        exit("ERROR: Fasta directory does not exist")
    app.config['fasta_dir'] = fasta_dir

    # Check that vizalign path exists
    if not os.path.isfile(vizalign):
        exit("ERROR: Invalid vizalign path")
    app.config['vizalign'] = vizalign

    app.register_blueprint(admin)
    Compress(app) # Compresses large responses (useful for alignment viewing)
    assets = Environment(app)
    js     = Bundle('js/d3.min.js',    
                    'js/colorbrewer.js',
                    'js/d3-tip.min.js',  
                    'js/panelutil.js',  
                    'js/scatterplot.js',
                    'js/jquery-2.1.3.min.js',
                    'js/jquery-ui.js',
                    'js/bootstrap.min.js',
                    'js/bootstrap-select.js',
                    'js/bubbleplot.js',
                    'js/barplot.js',
                    'js/nv.d3.js',
                    filters='jsmin', output='gen/packed.js')
    assets.register('str_validator_js_all', js)

    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///call_data.sqlite'
    app.config['ACTIVE_DB']               = HaploidCalls
    app.config['COMPARISONS']             = {"Father-son YSTRs" : HaploidCalls, "Marshfield" : DiploidCalls}
    app.config['ANALYSES']                = ["Bubble plot", "Bar plot"]
    app.config['MIN_BP_DIFF'] = -12
    app.config['MAX_BP_DIFF'] = 12
    db.init_app(app)
    with app.app_context():
        db.create_all()

    return app


if __name__== '__main__':
    parser = argparse.ArgumentParser(description='Web server to compare STR calls and their associated alignments')
    parser.add_argument("--bams",     type=str, required=True, help='Comma-separated list of bam files')
    parser.add_argument("--bais",     type=str, required=True, help='Comma-separated list of bam index files. Order must match that of the --bams arguments')
    parser.add_argument("--fasta",    type=str, required=True, help='Directory containing chromosome fasta files')
    parser.add_argument("--vizalign", type=str, required=True, help='Full path for vizalign executable')
    args = parser.parse_args()

    # Create app
    my_app = create_app(args.bams, args.bais, args.fasta, args.vizalign)

    # Run it
    my_app.run(host='0.0.0.0', port=6015, debug=True)
