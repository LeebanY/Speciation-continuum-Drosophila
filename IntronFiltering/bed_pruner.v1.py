#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage: bed_pruner.py -b <FILE> [-s <INT> -m <INT> -M <INT> -r <INT> -l <INT> -S <BOOL>] [-h]

  [Options]
    -b, --bed_f <FILE>	                       BED file (4th column must be shared ID!)
    -r, --seed <INT>                      	   Random seed [default: 12345]
    -s, --slop <INT>						   Slop X bases at beginning and end of each feature [default: 10]
    -m, --min_feature_length <INT>			   Minimum length of feature to be parsed [default: 200]
    -M, --max_feature_length <INT>			   Maximum length of feature to be parsed [default: 1000]
    -S, --skip_first_feature <BOOL>			   Skip first feature [default: True]
    -l, --region_length <INT>				   Length of regions to sample from within BED file per ID [default: 100]
    -h, --help                                 Show this message

"""

'''
[Goals]
- needs to write new BED
- plot length by ID
- needs to plot spans for resulting intervals

'''
from docopt import docopt
import pandas as pd
from tqdm import tqdm
import numpy as np
import sys

def intervals_to_sites(intervals):
    '''starts = np.array([0, 5, 8, 11, 15])
       ends   = np.array([2, 8, 9, 13, 18])
    clens : array([2, 5, 6, 8, 11])
    _sites: array([1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1])
    _sites: array([0, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1])
    _sites: array([1, 1, 4, 1, 1, 1,  3, 1, 3, 1, 1])
    sites : array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])'''
    starts, ends = intervals
    if starts is None or ends is None:
        return None
    clens = np.cumsum(ends - starts)
    if np.any(clens):
        sites = np.ones(clens[-1], dtype=np.int64)
        sites[0] = starts[0]
        sites[clens[:-1]] = starts[1:] - ends[:-1] + 1
        sites = sites.cumsum()
        return sites
    return None

def read_bed(bed_f):
	return pd.read_csv(bed_f, sep="\t", usecols=[0,1,2,5,9], names=['chrom', 'start', 'end', 'strand', 'name']).sort_values(
            ['chrom', 'start'], ascending=[True, True]).reset_index(drop=True)

def sample_regions(bed_f, region_length, max_feature_length, min_feature_length, skip_first_feature, slop, seed):
	bed_df = read_bed(bed_f)
	# slopping
	bed_df['start'] = bed_df['start'] + slop
	bed_df['end'] = bed_df['end'] - slop
	valid_idxs = bed_df['start'] < bed_df['end']
	bed_df = bed_df[valid_idxs]
	# length filters
	bed_df['length'] = bed_df['end'] - bed_df['start']
	bed_df = bed_df[(bed_df['length'] >= min_feature_length) & (bed_df['length'] <= max_feature_length)]

	groups = bed_df.groupby(by='name')
	rng = np.random.default_rng(seed)

	sampled_regions = []
	for name, df in tqdm(groups, total=len(groups), desc="[%] Sampling regions", ncols=100, disable=False):
		# skip first feature
		if skip_first_feature:
			df = df[1:]
		if not df.empty:
			#print(df)
			if not len(df['strand'].unique()) == 1:
				sys.exit("[X] The following gene has features on opposing strands.\n%s" % df)
			if not len(df['chrom'].unique()) == 1:
				sys.exit("[X] The following gene has features on different chroms.\n%s" % df)
			strand = df['strand'].unique()[0]
			chrom = df['chrom'].unique()[0]
			sites = intervals_to_sites((np.array(df['start']), np.array(df['end'])))
			if sites is not None and sites.shape[0] >= region_length:
				#print('sites', sites)
				random_start = 0 if sites.shape[0] == region_length else rng.integers(low=0, high=(sites.shape[0] - region_length))
				#print('sampled_sites = (', sites[random_start], sites[random_start + region_length], ")")
				sampled_sites = sites[random_start:random_start+region_length]
				#print('sampled_sites', sampled_sites, sampled_sites.shape)
				starts_ends = [(region[0], region[-1] + 1) for region in np.split(sampled_sites, np.where(np.diff(sampled_sites) > 1)[0]+1)]
				for start, end in starts_ends:
					sampled_regions.append([chrom, start, end, name, '.', strand])
	if len(sampled_regions):
		region_df = pd.DataFrame(sampled_regions)
		out_f = "%s.l_%s.m_%s.M_%s.s_%s.r_%s.pruned.bed" % (
			".".join(bed_f.split(".")[0:-1]),
			region_length,
			min_feature_length,
			max_feature_length,
			slop,
			seed
			)
		region_df.to_csv(out_f, sep="\t", index=False, header=False)


def sanitize_args(raw_args):
	args = {}
	args['--bed_f'] = raw_args['--bed_f']
	args['--max_feature_length'] = int(raw_args['--max_feature_length'])
	args['--min_feature_length'] = int(raw_args['--min_feature_length'])
	args['--region_length'] = int(raw_args['--region_length'])
	args['--seed'] = int(raw_args['--seed'])
	args['--skip_first_feature'] = bool(raw_args['--skip_first_feature'])
	args['--slop'] = int(raw_args['--slop'])
	return args

if __name__ == '__main__':
    __version__ = '0.1'
    args = sanitize_args(docopt(__doc__))
    print(args)
    sample_regions(
    	bed_f=args['--bed_f'],
    	region_length=args['--region_length'],
    	max_feature_length=args['--max_feature_length'],
    	min_feature_length=args['--min_feature_length'],
    	skip_first_feature=args['--skip_first_feature'],
    	slop=args['--slop'],
    	seed=args['--seed']
    	)
