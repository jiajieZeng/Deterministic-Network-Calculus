import copy
import csv
import math
import random
import time
import platform
import os
import psutil

import numpy as np
import networkx as nx
import pandas as pd
from scipy.optimize import linprog
from model import Model
from lpInstance import LPInstance
from candidateDesign import CandidateDesign
from annotatedDesign import AnnotatedDesign
from annotatedLP import AnnotatedLP
from predTopo import PredTopo

comp = []

# algo 1
def constructPartialNCPSTrajectory(C):

    T_C = [set() for i in range(len(list(C.graph.nodes)) + 5)]  # Set of time variables at input of si
    F_C = [set() for i in range(len(list(C.graph.nodes)) + 5)]  # set of flow variables at input of si
    E_C = set()  # Set of flow-specific entry time variables
    D_C = set()  # Set of flow-specific departure instances
    O = nx.DiGraph()  # Set of network-wide temporal orders
    flow_vis = [0 for i in range(len(list(C.graph.nodes)) + 5)]

    global comp

    # Resolving interdependency between flows.
    while True:
        cnt1 = 0
        for i in range(len(T_C)):
            cnt1 += len(T_C[i])

        # Introducing time variables.
        for j in range(len(C.paths)):  # 每一条流
            idx = len(C.paths[j]) - 1
            si = C.paths[j][idx]  # Tracing back from the sink node

            if flow_vis[j] != 1:
                v_name = 't' + str(C.time_variable_index)  # 递增id t1 t2 t3
                real_name = 'ARR(t^' + str(si) + '_' + str(j) + ')'

                C.time_variable_index += 1
                T_C[si].add(v_name)
                D_C.add(v_name)  # Departure instance of F_(j)

                C.departure_variable_map[v_name] = real_name
                C.departure_variable_rmap[real_name] = v_name
                C.time_variable_map[v_name] = real_name
                C.time_variable_rmap[real_name] = v_name

                flow_vis[j] = 1

            idx -= 1
            si = C.paths[j][idx]
            while idx != 0:
                for i in range(len(C.paths)):  # 每一条流
                    if si in C.paths[i]:  # 是否是经过si节点的流
                        sid = C.paths[i].index(si)
                        successor = C.paths[i][sid + 1]
                        for time_instance in T_C[successor]:
                            real_name = C.time_variable_map[time_instance]
                            x = 'ARR(t^' + str(successor) + '_' + str(i) + ')'
                            if real_name != x:
                                continue
                            if 'SC' in real_name:
                                continue

                            sct_u = 'SC(t^' + str(si) + '_' + str(i) + ')'
                            sct_v = 't' + str(C.time_variable_index)
                            flag = True
                            for t in T_C[si]:
                                if C.time_variable_map[t] == sct_u:
                                    flag = False
                                    break
                            if flag:
                                C.time_variable_index += 1
                                T_C[si].add(sct_v)
                                C.time_variable_map[sct_v] = sct_u
                                C.time_variable_rmap[sct_u] = sct_v

                            flag = True
                            art_u = 'ARR(t^' + str(si) + '_' + str(i) + ')'
                            art_v = 't' + str(C.time_variable_index)
                            for t in T_C[si]:
                                if C.time_variable_map[t] == art_u:
                                    flag = False
                                    break
                            if flag:
                                C.time_variable_index += 1
                                T_C[si].add(art_v)
                                C.time_variable_map[art_v] = art_u
                                C.time_variable_rmap[art_u] = art_v

                            sct_v = C.time_variable_rmap[sct_u]
                            art_v = C.time_variable_rmap[art_u]
                            O.add_edges_from([(sct_v, art_v), (art_v, time_instance), (sct_v, time_instance)])
                            if [sct_v, art_v] not in comp:
                                comp.append([sct_v, art_v])
                            if [art_v, time_instance] not in comp:
                                comp.append([art_v, time_instance])
                            if [sct_v, time_instance] not in comp:
                                comp.append([art_v, time_instance])
                            # print('({} <= {}), ({} <= {})'.format(sct_v, art_v, art_v, time_instance))
                idx -= 1
                si = C.paths[j][idx]

        cnt2 = 0
        for i in range(len(T_C)):
            cnt2 += len(T_C[i])

        if cnt1 == cnt2:
            break

    # Introducing entry-time variables.
    for j in range(len(C.paths)):
        si = C.paths[j][1]  # 相当于第一个节点，第0个是virtual source
        real_name = 'ARR(t^' + str(si) + '_' + str(j) + ')'
        v_name = C.time_variable_rmap[real_name]  # t1 t2 t3
        E_C.add(v_name)
        e_real_name = 'ARR(t_' + str(j) + ')'
        C.entry_variable_map[v_name] = e_real_name
        C.entry_variable_rmap[e_real_name] = v_name

    # Introducing flow variables
    for i in range(len(list(C.s_net))):
        si = C.s_net[i]
        tmp_set = T_C[si]
        flow_nums = []
        # 求出T^si和所有T^succ(si)的并集
        for j in range(len(C.paths)):
            if si in C.paths[j]:
                flow_nums.append(j)
                successor = C.paths[j][C.paths[j].index(si) + 1]
                tmp_set = tmp_set.union(T_C[successor])
        for v in tmp_set:  # F1 : F^(2)_(1)(t)    F2 : F^(2)_(1)(t)  123是全局单例
            # 对于每一条经过的流，每个时间变量都引入一个F
            for flow_num in flow_nums:
                f_v_name = 'F' + str(C.flow_variable_index)  # F1 F2 F3 F4 F5
                f_real_name = 'F^' + str(si) + '_' + str(flow_num) + '(' + str(v) + ')'  # F^(si)_(j)(v)
                flag = True
                for flow_variable in F_C[si]:
                    if C.flow_variable_map[flow_variable] == f_real_name:
                        flag = False
                if flag:
                    C.flow_variable_index += 1
                    F_C[si].add(f_v_name)
                    C.flow_variable_map[f_v_name] = f_real_name
                    C.flow_variable_rmap[f_real_name] = f_v_name

            for k in range(len(C.paths)):
                # 对于所有流的succ(si)都引入一个F
                if si in C.paths[k]:
                    k_successor = C.paths[k][C.paths[k].index(si) + 1]
                    fk_v_name = 'F' + str(C.flow_variable_index)
                    fk_real_name = 'F^' + str(k_successor) + '_' + str(k) + '(' + str(v) + ')'
                    flag = True
                    for flow_variable in F_C[k_successor]:
                        if C.flow_variable_map[flow_variable] == fk_real_name:
                            flag = False
                    if flag:
                        C.flow_variable_index += 1
                        F_C[k_successor].add(fk_v_name)
                        C.flow_variable_map[fk_v_name] = fk_real_name
                        C.flow_variable_rmap[fk_real_name] = fk_v_name

    # print(T_C, '\n', F_C, '\n', O.adj, '\n', E_C, '\n', D_C)
    C_prime = AnnotatedDesign(C, T_C, F_C, O, E_C, D_C)
    return C_prime

# algo 2
def constructCTOs(C_prime, filename):
    # edges = list(C_prime.O.edges())
    # for tup in edges:
    #     print(tup[0], ' ', tup[1])
    # return list(nx.all_topological_sorts(C_prime.O))
    # 从CSV文件读取DataFrame
    df = pd.read_csv(filename + 'predict_all_topo.csv', header=None)
    # 将DataFrame转换为二维列表
    topo = df.values.tolist()
    # print(topo)
    return topo

# algo 3
def completeNCPSTrajectory(C_prime, topo):
    B_T = set()  # Set of flow variable constraints
    time_constraints = set()
    for i in range(len(C_prime.C.s_net)):  # S_net
        si = C_prime.C.s_net[i]  # si
        # add service curve constraints
        for j in range(len(C_prime.C.paths)):  # 每一条流
            lhs = ''
            rhs = ''
            beta_rhs = ''
            if si not in C_prime.C.paths[j]:
                continue
            # si的后继节点
            successor = C_prime.C.paths[j][C_prime.C.paths[j].index(si) + 1]
            # 第j条流，经过si节点到达successor的时间
            lhs_real_name = 'ARR(t^' + str(successor) + '_' + str(j) + ')'
            lhs_v_name = C_prime.C.time_variable_rmap[lhs_real_name]
            # 左边求和
            for k in range(len(C_prime.C.paths)):  # 每一条流
                if si not in C_prime.C.paths[k]:
                    continue
                # 找到si这个节点在这条流的后继
                s_successor = C_prime.C.paths[k][C_prime.C.paths[k].index(si) + 1]
                lhs_flow_real_name = 'F^' + str(s_successor) + '_' + str(k) + '(' + lhs_v_name + ')'
                # 乘一个负号
                lhs = lhs + '-' + C_prime.C.flow_variable_rmap[lhs_flow_real_name]
            # 第j条流，到达si的后继时间t的service curve time
            rhs_real_name = 'SC(t^' + str(si) + '_' + str(j) + ')'
            rhs_v_name = C_prime.C.time_variable_rmap[rhs_real_name]
            for k in range(len(C_prime.C.paths)):
                if si not in C_prime.C.paths[k]:
                    continue
                # 采样
                rhs_flow_real_name = 'F^' + str(si) + '_' + str(k) + '(' + rhs_v_name + ')'
                rhs = rhs + '+' + C_prime.C.flow_variable_rmap[rhs_flow_real_name]
            rhs = rhs + '+' + '{:.10f}'.format(C_prime.C.SC[si][0]) + lhs_v_name + '-' + '{:.10f}'.format(C_prime.C.SC[si][0]) + rhs_v_name
            beta_rhs = beta_rhs + '+' + '{:.10f}'.format(C_prime.C.SC[si][0] * C_prime.C.SC[si][1])
            # >=
            # <=
            B_T.add(lhs + rhs + '<=' + beta_rhs)
            time_constraints.add('-' + lhs_v_name + '+' + rhs_v_name + '<=-' + '{:.10f}'.format(C_prime.C.SC[si][1]))

        # constraints on F^{[succ_j(s_i)]}_{(j)}(t_u)=F^{(s_i)}_{(j)}[FCFS(t_u)]
        for j in range(len(C_prime.C.paths)):
            if si not in C_prime.C.paths[j]:
                continue
            successor = C_prime.C.paths[j][C_prime.C.paths[j].index(si) + 1]
            lhs_real_name = 'ARR(t^' + str(successor) + '_' + str(j) + ')'
            lhs_v_name = C_prime.C.time_variable_rmap[lhs_real_name]
            rhs_real_name = 'ARR(t^' + str(si) + '_' + str(j) + ')'
            rhs_v_name = C_prime.C.time_variable_rmap[rhs_real_name]
            lhs_flow_real_name = 'F^' + str(successor) + '_' + str(j) + '(' + lhs_v_name + ')'
            lhs_flow_v_name = C_prime.C.flow_variable_rmap[lhs_flow_real_name]
            rhs_flow_real_name = 'F^' + str(si) + '_' + str(j) + '(' + rhs_v_name + ')'
            rhs_flow_v_name = C_prime.C.flow_variable_rmap[rhs_flow_real_name]
            B_T.add('+' + lhs_flow_v_name + '-' + rhs_flow_v_name + '=+0')

        # add causality constraints
        for j in range(len(C_prime.C.paths)):
            if si not in C_prime.C.paths[j]:
                continue
            successor = C_prime.C.paths[j][C_prime.C.paths[j].index(si) + 1]
            real_name = 'ARR(t^' + str(successor) + '_' + str(j) + ')'
            v_name = C_prime.C.time_variable_rmap[real_name]
            # si
            lhs_flow_real_name = 'F^' + str(si) + '_' + str(j) + '(' + v_name + ')'
            lhs_flow_v_name = C_prime.C.flow_variable_rmap[lhs_flow_real_name]
            # succj(si)
            rhs_flow_real_name = 'F^' + str(successor) + '_' + str(j) + '(' + v_name + ')'
            rhs_flow_v_name = C_prime.C.flow_variable_rmap[rhs_flow_real_name]
            B_T.add('-' + lhs_flow_v_name + '+' + rhs_flow_v_name + '<=+0')

            # time
            succ_start_u = 'ARR(t^' + str(si) + '_' + str(j) + ')'
            real = C_prime.C.time_variable_rmap[succ_start_u]
            lhs_flow_real_name = 'F^' + str(si) + '_' + str(j) + '(' + real + ')'
            lhs_flow_v_name = C_prime.C.flow_variable_rmap[lhs_flow_real_name]
            rhs_flow_real_name = 'F^' + str(successor) + '_' + str(j) + '(' + real + ')'
            rhs_flow_v_name = C_prime.C.flow_variable_rmap[rhs_flow_real_name]
            B_T.add('-' + lhs_flow_v_name + '+' + rhs_flow_v_name + '<=+0')

            # service cureve time
            sc_real_name = 'SC(t^' + str(si) + '_' + str(j) + ')'
            sc_v_name = C_prime.C.time_variable_rmap[sc_real_name]
            lhs_flow_real_name = 'F^' + str(si) + '_' + str(j) + '(' + sc_v_name + ')'
            lhs_flow_v_name = C_prime.C.flow_variable_rmap[lhs_flow_real_name]
            rhs_flow_real_name = 'F^' + str(successor) + '_' + str(j) + '(' + sc_v_name + ')'
            rhs_flow_v_name = C_prime.C.flow_variable_rmap[rhs_flow_real_name]
            B_T.add('-' + lhs_flow_v_name + '+' + rhs_flow_v_name + '<=+0')


    for i in range(len(list(C_prime.C.graph.nodes))):
        li = list(C_prime.C.graph.nodes)
        si = li[i]
        for t in C_prime.T_C[si]:
            for s in C_prime.T_C[si]:
                for j in range(len(C_prime.C.paths)):
                    # add monotonicity constraints
                    if si in C_prime.C.paths[j]:
                        t_id = topo.index(t)
                        s_id = topo.index(s)
                        lhs = C_prime.C.flow_variable_rmap['F^' + str(si) + '_' + str(j) + '(' + t + ')']
                        rhs = C_prime.C.flow_variable_rmap['F^' + str(si) + '_' + str(j) + '(' + s + ')']
                        if t_id > s_id:
                            B_T.add('-' + lhs + '+' + rhs + '<=+0')
                        elif t_id < s_id:
                            B_T.add('+' + lhs + '-' + rhs + '<=+0')
                        # add arrival curve constraints for F_(j)
                        if si == C_prime.C.paths[j][1]:
                            if t_id > s_id:
                                B_T.add('+' + lhs + '-' + rhs + '-' + '{:.10f}'.format(C_prime.C.flow_specific_AC[j][0]) + t
                                        + '+' + '{:.10f}'.format(C_prime.C.flow_specific_AC[j][0]) + s + '<=+' + '{:.10f}'.format(C_prime.C.flow_specific_AC[j][1]))
                            else:
                                B_T.add('+' + rhs + '-' + lhs + '-' + '{:.10f}'.format(C_prime.C.flow_specific_AC[j][0]) + s
                                        + '+' + '{:.10f}'.format(C_prime.C.flow_specific_AC[j][0]) + t + '<=+' + '{:.10f}'.format(C_prime.C.flow_specific_AC[j][1]))


        # # # add WCD constraints for topo
        # for j in range(len(C_prime.C.paths)):
        #     last_Pj = C_prime.C.paths[j][-1]
        #     last_real_name = 'ARR(t^' + str(last_Pj) + '_' + str(j) + ')'
        #     second_real_name = 'ARR(t' + '_' + str(j) + ')'
        #     last_v_name, second_v_name = '', ''
        #     for key, value in C_prime.C.departure_variable_map.items():
        #         if value == last_real_name:
        #             last_v_name = key
        #             break
        #     for key, value in C_prime.C.entry_variable_map.items():
        #         if value == second_real_name:
        #             second_v_name = key
        #             break
        #     time_constraints.add('+' + last_v_name + '-' + second_v_name + '<=+' + '{:.10f}'.format(C_prime.C.delay_constraints[j]))

    M_T = Model(C_prime.C, C_prime.T_C, C_prime.F_C, topo, B_T, C_prime.E_C, C_prime.D_C, time_constraints)
    return M_T

def saveAllTopo(topo):
    df = pd.DataFrame(topo)
    # 将DataFrame写入CSV文件
    df.to_csv('../all_topo.csv', index=False, header=False)  # 设置index和header为False以避免写入行和列的索引

# algo 4
def wholeNetworkWCDAnalysis(wholeCandidateSet, filename):
    LPs = []
    annotated_LPs = []
    global timeout_limit, timeout_start
    for value in wholeCandidateSet:
        C = value[0]
        C_id = value[1]
        time_constraints = set()
        C_prime = constructPartialNCPSTrajectory(C)  # algo1
        # print(C.time_variable_map)
        C_topo = constructCTOs(C_prime, filename)  # algo2
        # C_topo = []
        # saveAllTopo(C_topo)
        # constructNum = 0
        for topo in C_topo:
            objective = ''
            M_T = completeNCPSTrajectory(C_prime, topo)  # algo3
            time_constraints.clear()
            for j in range(len(C.delay_constraints)):
                # add flow-specific objective function for F_j
                departure_real_time = 'ARR(t^' + str(M_T.C.paths[j][-1]) + '_' + str(j) + ')'
                entry_real_time = 'ARR(t_' + str(j) + ')'
                departure_time, entry_time = '', ''
                for key, val in C.departure_variable_map.items():
                    if val == departure_real_time:
                        departure_time = key
                        break
                for key, val in C.entry_variable_map.items():
                    if val == entry_real_time:
                        entry_time = key
                        break
                objective = objective + '-' + departure_time + '+' + entry_time
            time_constraints.add(objective + '<=+0')

            # 在这里加入拓扑排序求得的时间大小关系
            for i in range(len(topo) - 1):
                time_constraints.add('+' + topo[i] + '-' + topo[i + 1] + '<=+0')

            lp = LPInstance(objective, M_T.B_T,
                                time_constraints.union(M_T.time_constraints),
                                topo, C, C_prime, C_id, C.time_variable_map,C.flow_variable_map,
                                None,None,None,None,None,None,None)  # add L_{T,C,j} into LPs
            objective_matrix, ne_matrix, ne_b, eq_matrix, eq_b = getMatrix(lp)
            lp.obj = objective_matrix
            lp.ne = ne_matrix
            lp.neb = ne_b
            lp.eq = eq_matrix
            lp.eqb = eq_b
            LPs.append(lp)
            # print(constructNum)
            # constructNum += 1
            # print(lp.time_map)
    # print(M_T.time_constraints)
    candidate_set_len = len(wholeCandidateSet)
    max_val = [0 for i in range(candidate_set_len)]
    all_li = [[] for i in range(candidate_set_len)]
    can_solve = [0 for i in range(candidate_set_len)]
    worst_case_cnt = [0 for i in range(candidate_set_len)]
    for L in LPs:
        timeout_start = time.time()
        # objective_matrix, ne_matrix, ne_b, eq_matrix, eq_b = getMatrix(L)
        x = linprog(L.obj, L.ne, L.neb, L.eq, L.eqb, method="highs")  # solve L with an LP solver to obtain L'
        # print(x)
        if x.success is True and (-1 * x.fun < 1):
            can_solve[L.id] += 1
            print(-1 * x.fun)
            L_prime = AnnotatedLP(-1 * x.fun, L.topo_list, L.C, L.C_prime, L.id,
                                    L.time_map, L.flow_map,
                                    find_new_edges(L.topo_list, L.C_prime.O),
                                    L.map, list(x.x), L.time_constraints, L.flow_constraints)
            annotated_LPs.append(L_prime)
            max_val[L.id] = max(max_val[L.id], L_prime.val)
    ans = []
    for L in annotated_LPs:
        if abs(L.val - max_val[L.C_id]) <= 1e-6 and worst_case_cnt[L.C_id] < 30:
            worst_case_cnt[L.C_id] += 1
            ans.append(L)
    for i in range(candidate_set_len):
        print('max:', i, max_val[i])
    return ans

def PredictwholeNetworkWCDAnalysis(wholeCandidateSet, pred_topo, filename):
    LPs = set()
    annotated_LPs = []
    # pred_val = [['case', 'val', 'topo']]
    pred_val = [['case', 'val']]
    wcd = [0.00009, 0.00009, 0.00009, 0.0001, 0.0003, 0.0004, 0.0003, 0.0003, 0.0003, 0.0003]
    for value in wholeCandidateSet:
        C = value[0]
        C_id = value[1]
        time_constraints = set()
        C_prime = constructPartialNCPSTrajectory(C)  # algo1
        C_T = pred_topo[C_id]
        print('------------------------------------------------------------------------------')
        print('case ', C_id, 'nums of prediction: ', len(C_T))
        for i in range(len(C_T)):
            topo = C_T[i]
            objective = ''
            M_T = completeNCPSTrajectory(C_prime, topo)  # algo3
            time_constraints.clear()
            for j in range(len(C.delay_constraints)):
                # add flow-specific objective function for F_j
                departure_real_time = 'ARR(t^' + str(M_T.C.paths[j][-1]) + '_' + str(j) + ')'
                entry_real_time = 'ARR(t_' + str(j) + ')'
                departure_time, entry_time = '', ''
                for key, value in C.departure_variable_map.items():
                    if value == departure_real_time:
                        departure_time = key
                        break
                for key, value in C.entry_variable_map.items():
                    if value == entry_real_time:
                        entry_time = key
                        break
                objective = objective + '-' + departure_time + '+' + entry_time
            time_constraints.add(objective + '<=+0')
            for i in range(len(topo) - 1):
                time_constraints.add('+' + topo[i] + '-' + topo[i + 1] + '<=+0')
            lp = LPInstance(objective, M_T.B_T, time_constraints.union(M_T.time_constraints), topo, C, C_prime, C_id, C.time_variable_map
                            , C.flow_variable_map,None,None,None,None,None,None,None)  # add L_{T,C,j} into LPs
            objective_matrix, ne_matrix, ne_b, eq_matrix, eq_b = getMatrix(lp)
            lp.obj = objective_matrix
            lp.ne = ne_matrix
            lp.neb = ne_b
            lp.eq = eq_matrix
            lp.eqb = eq_b
            x = linprog(objective_matrix, ne_matrix, ne_b, eq_matrix, eq_b, method='highs')
            if x.success is True:
                print(-1 * x.fun, topo)
                # pred_val.append([C_id, -1 * x.fun, topo])
                pred_val.append([C_id, -1 * x.fun])
            # else:
                # pred_val.append([C_id, 'None', topo])
                # pred_val.append([C_id, 'None'])
    writefile(filename + 'predict_result.csv', pred_val)

def isoper(c):
    flag = (c == '+') or (c == '-') or (c == '=') or (c == '<')
    return flag

def getMatrix(instance):
    id_map = dict()
    id_rmap = dict()
    id = 0
    for key, value in instance.time_map.items():
        id_map[key] = id
        id_rmap[id] = key
        id += 1
    for key, value in instance.flow_map.items():
        id_map[key] = id
        id_rmap[id] = key
        id += 1
    instance.map = copy.deepcopy(id_rmap)
    objective = np.zeros(id)
    i = 0
    while i < len(instance.objective):
        num = 0.
        t, oper, var = '', '+', ''
        if instance.objective[i] == '+':
            oper = '+'
        else:
            oper = '-'
        i += 1
        if instance.objective[i].isalpha():
            num = 1
        else:
            while i < len(instance.objective) and instance.objective[i].isdigit():
                t = t + instance.objective[i]
                i += 1
            num = float(t)
        while i < len(instance.objective) and (isoper(instance.objective[i]) is False):
            var += instance.objective[i]
            i += 1
        if oper == '-':
            num *= -1
        objective[id_map[var]] += num

    num_eq, num_ne = 0, 0
    for line in instance.flow_constraints:
        if '<=' in line:
            num_ne += 1
        else:
            num_eq += 1
    for line in instance.time_constraints:
        if '<=' in line:
            num_ne += 1
        else:
            num_eq += 1
    eq = np.zeros((num_eq, id))
    beq = np.zeros(num_eq)
    ne = np.zeros((num_ne, id))
    bne = np.zeros(num_ne)
    bound = np.zeros((id, 2))
    eq_row, ne_row = 0, 0
    all = instance.flow_constraints + instance.time_constraints
    for line in all:
        i = 0
        flag = True
        # 解析左边
        while flag:
            var = ''
            num = 0.
            oper = line[i]
            i += 1
            if line[i].isalpha():
                num = 1
            else:
                t = ''
                while (line[i].isalpha() is False) and (isoper(line[i]) is False):
                    t = t + line[i]
                    i += 1
                num = float(t)
            while isoper(line[i]) is False:
                var = var + line[i]
                i += 1

            if oper == '-':
                num *= -1
            if '<=' in line:
                ne[ne_row][id_map[var]] += num
            else:
                eq[eq_row][id_map[var]] += num
            if line[i] == '<' or line[i] == '=':
                flag = False
        if '<=' in line:
            i += 2
        else:
            i += 1
        num = 0.
        while i < len(line):
            oper, t = '+', ''
            if line[i] == '+' or line[i] == '-':
                oper = line[i]
                i += 1
            if line[i].isdigit():
                while i < len(line) and (isoper(line[i]) is False):
                    t += line[i]
                    i += 1
                tnum = float(t)
                if oper == '-':
                    tnum *= -1
                num += tnum
        if '<=' in line:
            bne[ne_row] += num
            ne_row += 1
        else:
            beq[eq_row] += num
            eq_row += 1
    return objective, ne, bne, eq, beq

def find_new_edges(topo_list, graph):
    new_edges = []
    for i in range(len(topo_list) - 1):
        lhs = topo_list[i]
        rhs = topo_list[i + 1]
        if rhs not in graph.adj[lhs]:
            new_edges.append((lhs, rhs))
    return new_edges

def getLabel(annotated, filename):
    lps = annotated
    nodes_csv = [['graph_id', 'node_id', 'label', 'feat']]
    nodes_map = [['map']]
    edges_csv = [['graph_id', 'src_id', 'dst_id']]
    graph_csv = [['graph_id', 'label']]
    topo_csv = [['graph_id', 'topo_list']]
    time_csv = [['graph_id', 'times']]
    linear_time_csv = [['case_id', 'times']]
    linear_topo_csv = [['case_id', 'topo']]
    linear_wcd_csv = [['case_id', 'wcd']]
    indent = 0
    for LP in lps:
        C = LP.C
        time_map = LP.time_map
        C_prime = LP.C_prime
        member_ship = []
        for i in range(len(C.paths)):
            code = ''
            for j in range(0, len(C.paths)):
                if j == i:
                    if len(code):
                        code += ','
                    code += '1'
                else:
                    if len(code):
                        code += ','
                    code += '0'
            member_ship.append(code)

        graph_id = indent
        f = indent * len(C_prime.O.nodes)
        topo_csv.append([graph_id, LP.topo_list])

        time_list = ''
        for idx in range(len(LP.digits)):
            name = LP.map[idx]
            if 't' not in name:
                continue
            if len(time_list):
                time_list += ','
            time_list += (name + '=' + '{:.10f}'.format(LP.digits[idx]))
        time_csv.append([graph_id, time_list])
        linear_time_csv.append([LP.C_id, time_list])
        indent += 1
        id_dif = len(lps)
        for key, value in time_map.items():
            id = key[1:]
            si, flow_num = 0, 0
            if 'SC' in value:
                si = int(value[5:value.index('_')])
                flow_num = int(value[value.index('_') + 1:value.index(')')])
            elif 'ARR' in value:
                si = int(value[6:value.index('_')])
                flow_num = int(value[value.index('_') + 1:value.index(')')])
            nodes_csv.append([graph_id, str(int(id) + f), str(LP.topo_list.index(key)),
                              '{}, {}, {:.18f}, {:.18f}, {:.18f}, {:.18f}'.format(graph_id / id_dif, member_ship[flow_num],
                                                                                  C.SC[si][0],
                                                                                  C.SC[si][1],
                                                                                  C.flow_specific_AC[flow_num][0],
                                                                                  C.flow_specific_AC[flow_num][1])])
            nodes_map.append([key])
        g = LP.C_prime.O
        for lhs, edge in g.adj.items():
            for rhs, value in edge.items():
                real_lhs = time_map[lhs]
                real_rhs = time_map[rhs]
                lhs_id = lhs[1:]
                rhs_id = rhs[1:]
                if 'SC' in real_lhs:
                    lhs_si = int(real_lhs[5:real_lhs.index('_')])
                    lhs_j = int(real_lhs[real_lhs.index('_') + 1: real_lhs.index(')')])
                elif 'ARR' in real_lhs:
                    lhs_si = int(real_lhs[6:real_lhs.index('_')])
                    lhs_j = int(real_lhs[real_lhs.index('_')+ 1: real_lhs.index(')')])
                if 'SC' in real_rhs:
                    rhs_si = int(real_rhs[5:real_rhs.index('_')])
                    rhs_j = int(real_rhs[real_rhs.index('_') + 1: real_rhs.index(')')])
                elif 'ARR' in real_rhs:
                    rhs_si = int(real_rhs[6:real_rhs.index('_')])
                    rhs_j = int(real_rhs[real_rhs.index('_') + 1: real_rhs.index(')')])
                edges_csv.append([graph_id, str(int(lhs_id) + f), str(int(rhs_id) + f)])

        graph_csv.append([graph_id, '{:.18f}'.format(LP.val)])
        linear_wcd_csv.append([LP.C_id, '{:.18f}'.format(LP.val)])
        linear_topo_csv.append([LP.C_id, LP.topo_list])
    # print('comp', comp)
    # print('nodes_csv', nodes_csv)
    # print('edges_csv', edges_csv)
    # print('graph_csv', graph_csv)
    # print('nodes_map', nodes_map)
    writefile(filename + 'nodes.csv', nodes_csv)
    writefile(filename + 'edges.csv', edges_csv)
    writefile(filename + 'graph.csv', graph_csv)
    writefile(filename + 'cons.csv', comp)
    writefile(filename + 'map.csv', nodes_map)
    writefile(filename + 'time.csv', time_csv)
    writefile(filename + 'topos.csv', topo_csv)
    writefile(filename + 'linear_topos.csv', linear_topo_csv)
    writefile(filename + 'linear_wcd.csv', linear_wcd_csv)
    writefile(filename + 'linear_time.csv', linear_time_csv)

def readfile(filename):
    # TODO:
    # add init
    # 这里需要
    # 读入物理图节点G，是列表形如[1, 2, 3, 4, 5, 6, 7, 8, 9]
    # 读入所有的流，每一条流是一个列表，形如 paths = [[0, 1, 4, 6, 8, 9], [0, 3, 5, 6, 9], [0, 2, 3, 6, 7, 9]]
    # 读入所有的S_net节点，形如s_net = [1, 2, 3, 4, 5, 6, 7, 8]
    # 读入所有的AC约束，形如(R,T)  flow_specific_AC = [(100, 10), (100, 10), (100, 10), (100, 10), (100, 10)]
    # 读入所有的SC约束，形如(p,delta)  SC = [(200, 20), (200, 20), (200, 20), (200, 20), (200, 20)]
    # 读入所有的D约束，形如   delay_constraints = [5]
    # 其他的暂时用不上，传入空列表占位
    # CSV input格式
    # 第一行，图的node
    # 第二行AC约束，一个p和一个delta 连续出现
    # 第三行SC约束，一个R和一个T 连续出现
    # 第四行D约束
    # 第五行流的数量
    # N条流
    data = []
    whole_candidate_set = list()
    with open(filename, 'r', newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            processed_row = []
            for item in row:
                try:
                    processed_row.append(float(item))
                except ValueError:
                    processed_row.append(item)
            data.append(processed_row)
    for ep in range(10):
        indent = ep * 10
        nodes = data[1 + indent]
        for i in range(len(nodes)):
            nodes[i] = int(nodes[i])
        s_net = data[2 + indent]
        for i in range(len(s_net)):
            s_net[i] = int(s_net[i])
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        flow_specific_AC = []
        for i in range(0, len(data[3 + indent]), 2):
            flow_specific_AC.append((data[3 + indent][i], data[3 + indent][i + 1]))
        SC = []
        for i in range(0, len(data[4 + indent]), 2):
            SC.append((data[4 + indent][i], data[4 + indent][i + 1]))
        wcd_constraints = []
        for i in range(len(data[5 + indent])):
            wcd_constraints.append(data[5 + indent][i])
        x = data[6 + indent][0]
        paths = []
        for i in range(int(x)):
            paths.append(data[7 + i + indent])
            for j in range(len(paths[i])):
                paths[i][j] = int(paths[i][j])
        C = CandidateDesign(G, flow_specific_AC, paths, SC, s_net, wcd_constraints)
        whole_candidate_set.append((C, ep))
    return whole_candidate_set

def writefile(filename, table):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(table[0])
        writer.writerows(table[1:])

def read_pred_file(filename):
    data = []
    for i in range(10):
        with open(filename + str(i) + '.csv', newline='') as file:
            reader = csv.reader(file)
            t_data = []
            for row in reader:
                t_data.append(row)
            data.append(t_data)
    return data

def run(filename, pred_filename, flag):
    whole_set = readfile(filename)
    filedir = filename[:filename.rfind('/')+1]
    if flag == 0:
        LPs = wholeNetworkWCDAnalysis(whole_set, filedir)
        getLabel(LPs, filedir)
    else:
        pred_topo = read_pred_file(pred_filename)
        PredictwholeNetworkWCDAnalysis(whole_set, pred_topo, filedir)

def show_info():
    # 计算消耗内存
    pid = os.getpid()
    # 模块名比较容易理解：获得当前进程的pid
    p = psutil.Process(pid)
    # 根据pid找到进程，进而找到占用的内存值
    info = p.memory_full_info()
    memory = info.uss / 1024 / 1024
    return memory

if __name__ == '__main__':
    start_memory = show_info()
    time_start = time.perf_counter()
    filename = '../test/202404/networkA3/case01.csv'
    pred_filename = '../test/202404/GNN-Result/networkA3-'
    run(filename, pred_filename, 1)
    time_end = time.perf_counter()
    time_sum = (time_end - time_start) * 1000
    print('time:{:.20f}'.format(time_sum))
    end_memory = show_info()
    print(f'memory: {end_memory - start_memory}MB')
