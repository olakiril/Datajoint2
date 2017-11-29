function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'obj', 'manolis_objects_old');
end
    obj = schemaObject;
end