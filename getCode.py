


def get_code(modName,tree, feature_name, n_features):
    left      = tree.tree_.children_left
    right     = tree.tree_.children_right
    threshold = tree.tree_.threshold
    feature_names=[feature_name+'[%i]'%i for i in range(n_features)]
    features  = [feature_names[i] for i in tree.tree_.feature]
    value = tree.tree_.value
    s1=0
    scode= 'void %s('%modName
    scode+= 'float *'+feature_name+','
    scode+='float *value){\n'
    def recurse(left, right, threshold, features, node, s1,scode):
        if (threshold[node] != -2):
            scode+=\
                "if ( " + features[node] + " <= " + str(threshold[node]) + " ) {"
            if left[node] != -1:
                s1,scode=recurse (left, right, threshold, features,left[node],s1,scode)
                scode+= "} else {\n"
                if right[node] != -1:
                    s1,scode=recurse (left, right, threshold, features,right[node],s1,scode)
                    scode+= "}\n"
        else:
                    #print "return " + str(value[node]/(tree.tree_.n_node_samples[node]+0.))
#            print "*p1=" + str(value[node][0][0]/(tree.tree_.n_node_samples[node]+0.))+';',
            scode+="*value="+str(value[node][0][0])+";\n"
#            print "// node_samples%=" + str(tree.tree_.n_node_samples[node]/(tree.tree_.n_node_samples[0]+0.))
            s1+=tree.tree_.n_node_samples[node]/(tree.tree_.n_node_samples[0]+0.)
            #print "*p2=" + str(value[node][0][1]/(tree.tree_.n_node_samples[node]+0.))+';'
        return s1,scode
    s1,scode=recurse(left, right, threshold, features, 0, s1,scode)
    #print 's1=',s1
    return scode+"}\n"
   

def get_codeEns(modName, feature_name,ntrees):
    scode= 'void %s('%modName
    scode+= 'float *'+feature_name+','
    scode+='float *value){\n'
    scode+="int i;\n"
    scode+="float vali;\n"
    scode+="*value=0;\n"
    for i in range(ntrees):
        fcall='%s%2.2i('%(modName,i)
        fcall+= feature_name+','
        fcall+='&vali);\n'
        scode+=fcall
        scode+="*value=*value+vali;\n"
        #scode+="printf(\"%g\\n\",vali);"
    scode+="*value=*value/%i;\n"%ntrees
    #scode+="printf(\"%g\\n\",*value);"
    scode+="}\n"
    return scode
