struct Saturation{D}
    dep::D
    is_saturated::Bool
end

function Base.isless(a::Saturation, b::Saturation)
    return isless((a.dep, a.is_saturated), (b.dep, b.is_saturated))
end

function MM.category_markershape(d::Saturation)
    return MM.category_markershape(d.dep)
end

function MM.category_markercolor(d::Saturation)
    return MM.category_markercolor(d.dep)
end

function MM.category_label(d::Saturation)
    label = MM.category_label(d.dep)
    if d.is_saturated
        label = "Saturated " * label
    end
    return label
end

function MM.category_markerstrokewidth(d::Saturation)
    if d.is_saturated
        return 0
    else
        return MM.category_markerstrokewidth(d.dep)
    end
end

function saturated_dependence(s::Iterator)
    dep = s.border.dependence
    return MM.BasisDependence(
        dep.basis,
        Saturation[
            Saturation(dep.dependence[i], is_saturated(s.matrix, mono))
            for (i, mono) in enumerate(dep.basis.monomials)
        ],
    )
end
