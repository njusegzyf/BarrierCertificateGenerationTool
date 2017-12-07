function res = isResNormsOk(resNorms)

    for norm = resNorms
        if norm >= 0.00001
            res = false;
            return;
        end
    end
    
    res = true;
end
