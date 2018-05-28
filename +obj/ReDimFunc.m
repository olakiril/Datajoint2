%{
# dimensionality reduction method
reduce_func       : varchar(256)     # function
%}
 
 
classdef ReDimFunc < dj.Lookup
     properties 
        contents = {
            '@(x,y) tsne(x'',[],y,50,30)'
            '@(x,y) lle(x,12,y)'''
            }
     end
end
