function q=func(y,Y,val)
q = sum((Y<y).*val)/sum(val);
end

