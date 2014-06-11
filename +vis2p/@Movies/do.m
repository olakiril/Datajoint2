function do( obj )
for key = fliplr(enumerate(obj))  % in reverse to show latest results first
    if isempty(MapJobs(key))
        jobkey = key;
        % reserve a job
        jobkey.status2 = 'reserved';
        populate( MapJobs(jobkey), 'quiet' );
        try
            % do the job
            populate(VisualMaps(key));

            % report job completed
            jobkey.status2 = 'completed';
            populate( MapJobs(jobkey), 'quiet', 'replace' );
        catch
            % report error
            err = lasterror;
            jobkey.status2 = 'error';
            jobkey.error_message = err.message;
            populate( MapJobs(jobkey), 'quiet', 'replace' );
        end
    end
end