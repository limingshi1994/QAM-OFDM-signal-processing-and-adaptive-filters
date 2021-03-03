function [BER] = ber(S, R)
   BER = 0;
   for i=1:length(S)
      if not(S(i) == R(i))
      BER = BER + 1;
      end
   end
   BER = BER*100/length(S);
   BER = sprintf("%.4f%%", BER);
end