function optcontract(ixs, xs, iy)
    size_dict = get_size_dict(ixs, xs)

    print("finding optimal tree ...")
    tree, cost = optimaltree(ixs, size_dict)
    @show tree, cost
    treecontract(tree, ixs, xs, iy)
end
