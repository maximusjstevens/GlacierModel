%%%% Take Al's Nisqually mass-balance reconstructions (from 11/2016) from  and reformat the data
%%%% so that it is easy to import into the model, i.e. csv w/ continuous
%%%% rows
%%%% Max Stevens, 12/2016

Bw=importdata('NisqBalance_winter.txt');
Bw2=reshape(Bw.data',1,[]);
w_years=1914:2015;
w_record=[w_years; Bw2(1:end-2)];

Bs=importdata('NisqBalance_summer.txt');
Bs2=reshape(Bs.data',1,[]);
s_years=1917:2015;
s_record=[s_years; Bs2(1:end-5)];

Ba=importdata('NisqBalance_annual.txt');
Ba2=reshape(Ba.data',1,[]);
a_years=1917:2015;
a_record=[a_years; Ba2(1:end-5)];

full_record=[a_years; Bw2(4:end-2); Bs2(1:end-5); Ba2(1:end-5)];

fnm='Nisq_MB_1917_2015.csv';
csvwrite(fnm,full_record)