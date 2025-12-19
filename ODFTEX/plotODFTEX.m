%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot ODFTEX textures with MTEX
%
%   Reference article:
%   
%   Ribe, N., Faccenda, M., VanderBeek, B.P. 
%   ODFTEX: A continuum model for texture evolution with 
%   dynamic recrystallization. Submitted to GJI, 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
%clf
%close all

path='add path here to file f1_ol.out';
addpath(path)

%Set colorscale limits
cmin = 0.0;
cstp = 0.25;
cmax = 6;

%Print
printmod = 0;

fontSize = 40;
setMTEXpref('FontSize',fontSize);
set(0,'DefaultAxesFontSize',fontSize);
set(0,'DefaultLegendFontSize',fontSize);

%Filename
fname0 = 'f_ol'; % plot olivine texture
%fname0 = 'f_opx'; % plot enstatite texture
fname1 = [fname0,'.out'];

% the crystal symmetry
cs = crystalSymmetry('222',[4.779 10.277 5.995]);
%cs = crystalSymmetry(l'orthorhombic');
ss = specimenSymmetry('triclinic');

% Compute ODF
odf = loadODF_generic(fname1,'CS',cs,...
    'ColumnNames',{'Euler1','Euler2','Euler3','weights'},...
    'kernel',SO3DeLaValleePoussinKernel('halfwidth',10*degree)); %'interp'

% Plot pole figures
plotx2east
h = Miller({1,0,0},{0,1,0},{0,0,1},cs,'uvw');
figure(12)
plotPDF(odf,h,'colorRange',[cmin,cmax],'contourf',cmin:cstp:cmax,'FontSize',20);
%plotPDF(odf,h,'3d','FontSize',20)
mtexColorbar

if printmod
    figname = [fname0,'_max',num2str(cmax),'.png'];
    exportgraphics(gcf, figname,'Resolution',300)
    movefile(figname,path)
end