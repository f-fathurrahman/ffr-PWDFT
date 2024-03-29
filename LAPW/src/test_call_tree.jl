mutable struct Node{T}
    data::T
    parent::Union{Nothing,Node{T}}
    children::Union{Nothing,Vector{Node{T}}}
    #
    function Node{T}(data::T; parent=nothing, children=nothing) where T
        new{T}(data, parent, children)
    end
end


function AbstractTrees.children(node::Node)
    if isnothing(node.children)
        return ()
    else
        return node.children
    end
end

function AbstractTrees.nodevalue(node::Node)
    return node.data
end


function get_called_subs(filename::String)
    f = open(filename)
    lines = readlines(f)
    close(f)

    called_subroutines = Vector{String}()
    pattern = r"\W*[Cc][Aa][Ll][Ll]\W*([A-Za-z0-0_]*)."

    for l in lines
        res = match(pattern, l)
        if !isnothing(res)
            push!(called_subroutines, String(res.captures[1]))
        end
    end
    return unique(called_subroutines)
end

my_tree = Node{String}("gndstate")
children = Vector{Node{String}}()
for s in called_subroutines
    push!(children, Node{String}(s, parent=my_tree))
end

