function C = int_coscos(m,n,a,b)
% ∫_a^b cos(mπx) cos(nπx) dx

C = 0.5 * ( int_cos_kpi(m-n,a,b) + int_cos_kpi(m+n,a,b) );
end
