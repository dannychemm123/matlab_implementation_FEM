function I = int_cos_kpi(k,a,b)
% ∫_a^b cos(kπx) dx

if k == 0
    I = b - a;
else
    I = (sin(k*pi*b) - sin(k*pi*a)) / (k*pi);
end
end
