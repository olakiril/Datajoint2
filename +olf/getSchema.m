function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'olf', 'manolis_olfaction');
end
obj = schemaObject;
end
