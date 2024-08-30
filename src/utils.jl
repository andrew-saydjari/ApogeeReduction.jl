ENV["SLACK_CHANNEL"]="C07KQ7BJY5P"

function initalize_git(git_dir)
    git_commit = LibGit2.head(git_dir)
    git_repo = LibGit2.GitRepo(git_dir)
    git_head = LibGit2.head(git_repo)
    git_branch = LibGit2.shortname(git_head)
    println("Running on branch: $git_branch, commit: $git_commit"); flush(stdout)
    return git_branch, git_commit
end

function isnanorzero(x)
    return isnan(x) | iszero(x)
end

nanzeromean(x) = if all(isnanorzero,x)
    NaN
else
    mean(filter(!isnanorzero,x))
end
nanzeromean(x,y) = mapslices(nanzeromean,x,dims=y)

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

nanzerosum(x) = if all(isnanorzero,x)
    NaN
else
    sum(filter(!isnanorzero,x))
end
nanzerosum(x,y) = mapslices(nanzerosum,x,dims=y)

nanzeromedian(x) = if all(isnanorzero,x)
    NaN
else
    median(filter(!isnanorzero,x))
end
nanzeromedian(x,y) = mapslices(nanzeromedian,x,dims=y)

nanzeroiqr(x) = if all(isnanorzero,x)
    NaN
else
    iqr(filter(!isnanorzero,x))/1.34896
end
nanzeroiqr(x,y) = mapslices(nanzeroiqr,x,dims=y)