function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'oldmap', 'manolis_map');
end
obj = schemaObject;
end
