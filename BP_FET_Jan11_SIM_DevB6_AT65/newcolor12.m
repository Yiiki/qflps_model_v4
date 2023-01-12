function colorPck=newcolor12()

stringA=["#4e3fee"
    "#1f78b4"
    "#b2df8a"
    "#33a02c"
    "#fb9a99"
    "#e31a1c"
    "#fdbf6f"
    "#ff7f00"
    "#cab2d6"
    "#6a3d9a"
    "#f8309f"
    "#b12828"];
charA=convertStringsToChars(stringA);
numA=zeros(size(charA,1),3);

for i=1:size(charA,1)
numA(i,:)=hex2rgb(charA{i});
end

color12=numA;


stringA=stringA(1:end-1);
charA=convertStringsToChars(stringA);
numA=zeros(size(charA,1),3);

for i=1:size(charA,1)
numA(i,:)=hex2rgb(charA{i});
end

color11=numA;
color11 = [0 0.4470 0.7410
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330
             0.6350 0.0780 0.1840
             0.25 0.80 0.54
             0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80];

color3=[240,128,128
    205,92,92
    220,20,60
    ]./256;

colorPck.c3=color3;
colorPck.c11=color11;
colorPck.c12=color12;

end