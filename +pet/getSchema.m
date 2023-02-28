function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'pet', 'petros_test');
end
obj = schemaObject;
end