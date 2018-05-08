function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'sim', 'manolis_simfilters');
end
    obj = schemaObject;
end