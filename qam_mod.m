function [QAM] = qam_mod(Sequence, M)
   QAM = qammod(Sequence, M, 'InputType', 'bit', 'UnitAveragePower',true);
end