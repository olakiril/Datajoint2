function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'mei', 'neurostatic_mei_closed_loop');
end
obj = schemaObject;
end
