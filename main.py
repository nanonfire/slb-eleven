# coding=utf-8
# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

# Read PCR data into a pandas DataFrame. You want a data file where each
# row corresponds to a separate well, with columns for the sample name,
# target name, and Cq value. NTC wells should have the sample name set to
# a value like 'NTC'.

import eleven

df = pd.read_csv('TestData.csv')

# If your Sample, Target, and Cq columns are called other things, they
# should be renamed to Sample, Target, and Cq.
df = df.rename(columns={'Gene': 'Target', 'Ct': 'Cq'})
df = df[df['Rat'].map(len) < 1]

# Drop the wells that are too close to the NTC for that target.
censored = eleven.censor_background(df)

# Rank your candidate reference genes.
ranked = eleven.rank_targets(censored, ['CypA', 'Rpl32', 'Hprt1'], 'Control')

# Normalize your data by your most stable genes and compute normalization
# factors (NFs).
nf = eleven.calculate_nf(censored, ranked.ix['Target', 0:3], 'Control')

# Now, normalize all of your expression data.
censored['RelExp'] = eleven.expression_nf(censored, nf, 'Control')
