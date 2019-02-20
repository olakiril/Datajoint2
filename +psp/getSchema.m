function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'psp', 'manolis_federico_events');
end
obj = schemaObject;
end
