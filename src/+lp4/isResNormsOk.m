function res = isResNormsOk(resNorms)

    import lp4.Lp4Config
    for norm = resNorms
        if norm >= Lp4Config.RES_NORM_THRESHOLD
            res = false;
            return;
        end
    end
    
    res = true;
end
