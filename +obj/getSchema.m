function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'obj', 'manolis_objects');
end
    obj = schemaObject;
end