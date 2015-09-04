function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    vis2p.getSchema
    schemaObject = dj.Schema(dj.conn, 'beh', 'manolis_behavior');
end
obj = schemaObject;
end
