clear;close all;

tt=struct;

tt.y2002.elevation=[nan	nan	2104	1887	nan	1720];
tt.y2002.winter=[nan	nan	3.16	2.72	nan	2.64];
tt.y2002.summer=[nan	nan	-5.01	-4.98	nan	-3.50];
tt.y2002.annual=[nan	nan	-1.85	-2.26	nan	-0.86];

tt.y2003.elevation=[nan	3025	2250	1875	1875	1750];
tt.y2003.winter=[nan	1.80	3.17	2.70	2.54	2.72];
tt.y2003.summer=[nan	-3.99	-6.45	-6.48	-5.43	-9.00];
tt.y2003.annual=[nan	-2.19	-3.28	-3.78	-2.89	-7.85];

tt.y2004.elevation=[nan	3018	2175	1890	1865	1775];
tt.y2004.winter=[nan	2.27	3.83	2.50	2.10	1.70];
tt.y2004.summer=[nan	-2.55	-4.97	-5.47	-3.18	nan];
tt.y2004.annual=[nan	-0.27	-1.14	-2.98	-1.07	nan];

tt.y2005.elevation=[nan	3020	2180	1875	1875	1750 ];
tt.y2005.winter=[nan	1.50	1.99	1.39	1.00	0.50 ];
tt.y2005.summer=[nan	-2.63	-5.15	-6.99	-2.02	-5.60 ];
tt.y2005.annual=[nan	-1.13	-3.16	-5.60	-1.01	-5.09 ];

tt.y2006.elevation=[3382	2960	2180	1875	1878	1774 ];
tt.y2006.winter=[2.40	2.67	3.66	3.01	3.71	2.99 ];
tt.y2006.summer=[-1.98	-2.45	-5.77	-7.82	-5.65	-6.90 ];
tt.y2006.annual=[0.41	0.22	-2.11	-4.81	-1.94	-3.92 ];

tt.y2007.elevation=[3382	2960	2175	1890	1870	1778 ];
tt.y2007.winter=[5.48	2.34	4.20	2.30	2.23	1.90 ];
tt.y2007.summer=[-2.80	-3.40	-5.46	-7.31	-4.36	-5.83 ];
tt.y2007.annual=[2.68	-1.06	-1.26	-5.01	-2.13	-3.92 ];

tt.y2008.elevation=[3382	2960	2175	1890	1870	1778 ];
tt.y2008.winter=[3.97	2.28	3.78	3.68	3.14	2.90 ];
tt.y2008.summer=[-1.80	-3.20	-5.16	-5.77	-4.97	-6.12 ];
tt.y2008.annual=[2.17	-0.92	-1.38	-2.09	-1.83	-3.22 ];

tt.y2009.elevation=[3382	2960	2175	1890	1870	1778 ];
tt.y2009.winter=[1.99	2.25	2.61	2.14	1.47	1.51 ];
tt.y2009.summer=[-2.33	-3.78	-6.24	-6.73	-2.53	-4.77 ];
tt.y2009.annual=[-0.34	-1.53	-3.64	-4.58	-1.06	-3.26 ];

tt.y2010.elevation=[3382	2960	2175	1890	1870	1778 ];
tt.y2010.winter=[3.40	3.31	3.13	2.45	2.50	2.24 ];
tt.y2010.summer=[-1.57	-1.75	-3.17	-4.34	-3.55	-4.82 ];
tt.y2010.annual=[1.83	1.56	-0.05	-1.89	-1.05	-2.57 ];

tt.y2011.elevation=[3382	2960	2175	1890	1870	1778 ];
tt.y2011.winter=[3.48	3.04	3.55	3.76	4.01	2.64 ];
tt.y2011.summer=[-1.95	-2.66	-4.71	-4.57	-3.64	-4.54 ];
tt.y2011.annual=[1.53	0.38	-1.16	-0.80	0.37	-1.89 ];

tt.y2012.elevation=[3382	2960	2175	1890	1870	1778 ];
tt.y2012.winter=[2.64	2.92	3.43	3.31	2.88	2.14 ];
tt.y2012.summer=[-2.30	-2.86	-4.52	-6.19	-2.96	-3.65 ];
tt.y2012.annual=[0.34	0.06	-1.10	-2.89	-0.07	-1.51 ];

years=2002:2012;




figure(1);
hold on;
grid on;
box on;
for ii=1:length(years)
    
    year=years(ii);
    yr=num2str(year);
    xx=strcat('y',yr);
    plot(tt.(xx).winter,tt.(xx).elevation)
end

figure(2);
hold on;
grid on;
box on;
for ii=1:length(years)
    
    year=years(ii);
    yr=num2str(year);
    xx=strcat('y',yr);
    plot(tt.(xx).summer,tt.(xx).elevation)
end

figure(3);
hold on;
grid on;
box on;
for ii=1:length(years)
    
    year=years(ii);
    yr=num2str(year);
    xx=strcat('y',yr);
    plot(tt.(xx).annual,tt.(xx).elevation)
end










