function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    vis2p.getSchema
    schemaObject = dj.Schema(dj.conn, 'oldbeh', 'manolis_behavior');
end
obj = schemaObject;
end
