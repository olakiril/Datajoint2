function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'simf', 'manolis_simfilters');
end
    obj = schemaObject;
end