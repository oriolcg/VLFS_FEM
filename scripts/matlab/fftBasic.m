function [fHalf,fAmp,fFull,fPhi] = fftBasic(fnc,sampf)
    
    %% fnc is 1d only
    
    sz=max(size(fnc));
    if(mod(sz,2)~=0)
        sz=sz-1;
        fnc=fnc(1:end-1);
    end
    
    fprintf('Sample Len = %d \n',sz);
    fprintf('Least count Hz = %f \n',sampf/sz);
    fprintf('Max Freq (Half band) Hz = %f \n \n',sampf/2);
    
    fAmp=fft(fnc);
    fPhi=unwrap(angle(fAmp));
    fAmp=abs(fAmp/sz);
    fAmp=fAmp(1:sz/2+1);
    fAmp(2:end-1)=2*fAmp(2:end-1);
    fFull=sampf*[0:sz-1]/sz;
    fHalf=sampf*[0:(sz/2)]/sz;
    %fHalf=fFull;

end