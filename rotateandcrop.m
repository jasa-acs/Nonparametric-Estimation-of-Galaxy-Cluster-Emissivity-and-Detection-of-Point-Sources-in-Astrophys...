function [F,E,O,Fp,Ep,Op]=rotateandcrop(F1,E1,O1,newn,center,angle)
    F=rotateandcrop1(F1,center,angle,newn);
    E=rotateandcrop1(E1,center,angle,newn);
    O=rotateandcrop1(O1,center,angle,newn);
    imagesc(log(F))
    Fp=fliplr(F);
    Ep=fliplr(E);
    Op=fliplr(O);
end

function A=rotateandcrop1(A,center,angle,newn)
    A=imrotate(A,angle,'nearest');
    if ~isempty(center)
        A=A((center(1)-newn/2+1):(center(1)+newn/2),(center(2)-newn/2+1):(center(2)+newn/2));
    end
    
end

function pepe=pepepotamo(s)
    ind=5
    realdata(ind).F=F;
    realdata(ind).E=E;
    realdata(ind).O=O;
    realdata(ind).Fp=Fp;
    realdata(ind).Ep=Ep;
    realdata(ind).Op=Op;
    realdata(ind).telescope='xmm';
    realdata(ind).name='perseus';
end