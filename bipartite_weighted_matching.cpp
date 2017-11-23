#include<iostream>
#include<type_traits>
#include<functional>

//find a min/max weight cardinality k matching of a dense bipartite graph with n1 left vertices and n2 right vertices
//set PR std::less if you want a minimum weight matching, std::greater if you want a maximum one
//c is an n1*n2 cost matrix
//mu1 (mu2) is an array whose i-th element is the partner of i-th left (right) vertex (-1 if unmatched)
template<template<typename>class PR = std::less, class mat, class vec1, class vec2>
void bipartite_weighted_matching(const mat &c, const int n1, const int n2, int k, vec1 &mu1, vec2 &mu2) {
	typedef typename remove_const<typename remove_reference<decltype(c[0][0])>::type>::type number;
	number *pi = new number[n2];
	int *pred = new int[n2];
	PR<number> cmp;
	auto dist = [&](int v) {
		int u = pred[v], w = mu1[u];
		return w == -1 ? c[u][v] : pi[w] - c[u][w] + c[u][v];
	};
	for (int i = 0;i < n1;i++)mu1[i] = -1;
	for (int j = 0;j < n2;j++)mu2[j] = -1;
	for (int i = 0;i < n1;i++)for (int j = 0;j < n2;j++)if (i == 0 || cmp(c[i][j], pi[j]))pi[j] = c[i][j], pred[j] = i - n1;
	for (int v = -1;;) {
		for (int j = 0;j < n2;j++)if (mu2[j] == -1 && (v == -1 || cmp(pi[j], pi[v])))v = j;
		while (v != -1) {
			int u = pred[v] + n1, w = mu1[u];
			mu2[v] = u;
			mu1[u] = v;
			v = w;
		}
		if (--k == 0)return delete[]pi, delete[]pred;
		for (int i = 0;i < n1;i++) {
			if (mu1[i] != -1)continue;
			for (int j = 0;j < n2;j++)if (pred[j] < 0 || cmp(c[i][j], c[pred[j]][j]))pred[j] = i;
		}
		for (;;v = -1) {
			for (int j = 0;j < n2;j++)if (pred[j] >= 0 && (v == -1 || cmp(dist(j) - pi[j], dist(v) - pi[v])))v = j;
			if (v == -1)break;
			pi[v] = dist(v);
			pred[v] -= n1;
			int u = mu2[v];
			if (u != -1)for (int j = 0;j < n2;j++)if (pred[j] >= 0 && cmp(pi[v] - c[u][v] + c[u][j], dist(j)))pred[j] = u;
		}
	}
}
