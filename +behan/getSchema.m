function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'behan', 'manolis_beh_analysis');
end
obj = schemaObject;
end
