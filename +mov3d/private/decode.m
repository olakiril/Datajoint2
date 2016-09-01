function mi = decode(data,dec_method)
switch dec_method
    case 'nnclassRaw'
        mi = nnclassRaw(data);
    case 'nnclassRawSV'
        mi = nnclassRawSV(data);
    case 'nnclass'
        mi = nnclass(data);
    case 'nnclassSV'
        mi = nnclassSV(data);
end
end