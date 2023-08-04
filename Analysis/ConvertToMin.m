function minObjFcn = ConvertToMin(x, objfcn)
    fcn = objfcn(x);
    if fcn >= 0
        minObjFcn = 1/(1+fcn);
    else
        minObjFcn = 1 + abs(fcn);
    end
  