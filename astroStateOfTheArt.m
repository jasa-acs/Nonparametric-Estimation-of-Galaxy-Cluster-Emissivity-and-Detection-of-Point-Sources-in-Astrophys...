function sol=astroStateOfTheArt(F,E,O, angle)
%set angle p.ex  [0 180], set [0 360] for full inversion
%compile first deprojection.C

M=size(F,1);
fmt=[repmat('%f ',1,M-1) '%f\n'];

fileID = fopen('Fcursim','w');
fprintf(fileID,fmt,F);
fclose(fileID);

fileID = fopen('Ecursim','w');
fprintf(fileID,fmt,E);
fclose(fileID);

fileID = fopen('Ocursim','w');
fprintf(fileID,fmt,O);
fclose(fileID);
Mbin=M/2;

rundeprojection=['deprojection Fcursim Ecursim Ocursim ',num2str(Mbin),' ',num2str(Mbin),' ',num2str(angle(1)),' ',num2str(angle(2)),' outcur']; %WARNING: SECTOR DEPROJECTION'

system(rundeprojection);
result = dlmread('outcur.txt','',1,0);
xbin=result(:,1);
fbin=result(:,4);
errorbin=result(:,3);
sol.xbin=xbin;
sol.fbin=fbin;
sol.errorbin=errorbin;
end

