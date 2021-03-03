function [out_aligned] = alignIO(out,pulse,IRlength)
    fs = 16e3;
    cross_corr = xcorr(out,pulse);
    cross_corr = cross_corr(length(out):end);
    [max_cross_corr,index] = max(cross_corr);
    
 
    out_aligned = out(index:end);
    cut = length(pulse)+IRlength-200;
    out_aligned = out_aligned(cut:end);
    
    figure(5)
    plot(out_aligned)
    title('aligned')
    figure(6)
    plot(out)
    title('out')
    
end

