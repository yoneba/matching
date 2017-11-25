#include<iostream>
#include<type_traits>
#include<functional>
#include<utility>

//find a min/max weight cardinality lambda matching of a complete bipartite graph with m left vertices and n right vertices
//set PR std::less if you want a minimum weight matching, std::greater if you want a maximum one
//c is an m*n cost matrix
//mu (nu) is a resulting array whose i-th element is the partner of i-th left (right) vertex (-1 if unmatched)
template<template<typename>class PR = std::less, class mat, class vec1, class vec2>
void bipartite_weighted_matching(const mat &c, const int m, const int n, const int lambda, vec1 &mu, vec2 &nu) {
	typedef typename remove_const<typename remove_reference<decltype(c[0][0])>::type>::type number;
	number *pi = new number[n], *dist = new number[n];
	int *pred = new int[n];
	PR<number> cmp;
	for (int i = 0;i < m;i++)mu[i] = -1;
	for (int j = 0;j < n;j++)nu[j] = -1;
	for (int k = 0;k < lambda;k++) {
		for (int j = 0;j < n;j++)pred[j] = -1;
		for (int i = 0;i < m;i++)if (mu[i] == -1)for (int j = 0;j < n;j++)if (pred[j] == -1 || cmp(c[i][j], dist[j]))dist[j] = c[i][j], pred[j] = i;
		for (int h = 0;h < k;h++) {
			int v = -1;
			for (int j = 0;j < n;j++)if (nu[j] != -1 && pred[j] >= 0 && (v == -1 || cmp(dist[j] - pi[j], dist[v] - pi[v])))v = j;
			pi[v] = dist[v], pred[v] -= m;
			int u = nu[v];
			for (int j = 0;j < n;j++)if (pred[j] >= 0 && cmp(pi[v] - c[u][v] + c[u][j], dist[j]))dist[j] = pi[v] - c[u][v] + c[u][j], pred[j] = u;
		}
		int v = -1;
		for (int j = 0;j < n;j++)if (nu[j] == -1 && (v == -1 || cmp(dist[j], dist[v])))v = j;
		if (v != -1)pi[v] = dist[v], pred[v] -= m;
		while (v != -1) {
			int u = pred[v] + m;
			nu[v] = u;
			std::swap(mu[u], v);
		}
	}
	delete[]pi, delete[]dist, delete[]pred;
}
