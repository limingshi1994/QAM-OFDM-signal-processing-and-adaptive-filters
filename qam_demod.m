function [Sequence] = qam_demod(QAM, M)
   Sequence = qamdemod(QAM,M,'OutputType','bit','UnitAveragePower',true);
end