unByteCode <- function(fun)
    {
        FUN <- eval(parse(text=deparse(fun)))
        environment(FUN) <- environment(fun)
        FUN
    }

assignEdgewise <- function(name, env, value)
    {
        unlockBinding(name, env=env)
        assign( name, envir=env, value=value)
        lockBinding(name, env=env)
        invisible(value)
    }

unByteCodeAssign <- function(fun)
    {
        name <- gsub('^.*::+','', deparse(substitute(fun)))
        FUN <- unByteCode(fun)
        retval <- assignEdgewise(name=name,
                                 env=environment(FUN),
                                 value=FUN
                                 )
        invisible(retval)
    }
