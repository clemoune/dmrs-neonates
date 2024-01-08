function ADCsphere = compute_ADC_from_Sphere(p, protocol)

S = SynthMeasSphere(p, protocol);
bvals = GetBvalues(protocol)./1E9;
ADCsphere = [(log( S(1) ) - log( S(2) ))/(bvals(2) - bvals(1)), ...
             (log( S(3) ) - log( S(4) ))/(bvals(4) - bvals(3)), ...
             (log( S(5) ) - log( S(6) ))/(bvals(6) - bvals(5)), ...
             (log( S(7) ) - log( S(8) ))/(bvals(8) - bvals(7))];

end
    
    
    