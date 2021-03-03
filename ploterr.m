function ploterr(err, FC, Lt)
    figure(100)
    for i = 1:(FC-Lt)
        err_max(i) = max(err(:,i));
    end
    plot(1:1:FC-Lt, err_max)
end
    