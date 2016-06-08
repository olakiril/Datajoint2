function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'mov3d', 'manolis_movies3d');
end
    obj = schemaObject;
end