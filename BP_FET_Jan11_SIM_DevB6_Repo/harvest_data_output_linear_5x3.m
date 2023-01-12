% harvest_data.m [followed]

color11=[
    56, 145, 204
    219, 91, 36
    237, 178, 35
    134, 60, 149
    124, 175, 55
    79, 191, 238
    162, 20, 47
    0, 114, 189
    199, 64, 33
    243, 200, 98
    81, 158, 209
    ]./256;

color11_cmp=1-color11;

color22=[color11;color11_cmp];

set(groot,'defaultAxesColorOrder',color11,...
      'defaultAxesLineStyleOrder','-|--|:')

sub_tag='a':'z';

idxSet=[
    3,1
    3,2
    3,4
    5,3
    6,5
    ];

fetnum=size(idxSet,1);

output_vds=linspace(0,3,101)';
output_vgs=linspace(-3,3,11);

pos_left=1;% cm
pos_botm=2;% cm
pos_widt=18.6;% cm
pos_heig=24.5;% cm

pos_sub_left=0.1;
pos_sub_botm=0.855;
pos_sub_spac_horizont=0.33;
pos_sub_spac_vertical=0.20;

sub_cols=3;
sub_rows=5;

h1=figure;
h1.Units='centimeters';
h1.Position=[pos_left pos_botm pos_widt pos_heig];

for i=1:fetnum

cols=mod(i-1,sub_cols)+1;
rows=floor((i-1)/sub_cols)+1;

sub_pos1=pos_sub_left+(cols-1)*pos_sub_spac_horizont;
sub_pos2=pos_sub_botm-(rows-1)*pos_sub_spac_vertical;
sub_pos3=0.2;
sub_pos4=0.115;

sub_pos_vec=[sub_pos1,sub_pos2,sub_pos3,sub_pos4];

fprintf('%d\n',i)

idsbco=ids_output_dat{idxSet(i,1),idxSet(i,2)};
output_ids=ids_output_sim{idxSet(i,1),idxSet(i,2)};

sh1=subplot('Position',sub_pos_vec);
plot(output_vds,idsbco.*1e5,'o',output_vds,output_ids.*1e5,'-','MarkerSize', 4)
stag=sprintf('(%s) Dev#%d-%d',sub_tag(i),idxSet(i,:));
title(stag)
xlim([0 3])
ylim([-0.05 1.0])
yticks([0 0.2 0.4 0.6 0.8 1.0])
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'})

xlabel('{\it V_{ds}} (V)')
ylabel('{\it I_{ds}} (10 uA)')

sh1.FontName='Arial';
sh1.FontSize=10;
sh1.LineWidth=1;

pause(0.1)

end

set(groot,'defaultAxesLineStyleOrder','remove')

set(groot,'defaultAxesColorOrder','remove')
