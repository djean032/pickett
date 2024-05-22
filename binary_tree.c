/*

Creator: Dairen Jean 22 May 2024

Organization: United States Military Academy, West Point, NY, 10996

Description:
    Binary Search Tree Implementation for SPFIT
    Used to compare exp. frequencies and calc. frequencies to determine nearest
    neighbor and assign appropriate quantum numbers. Used primarily for coupled
    vibrational states that exhibit state mixing.

*/

#include "calpgm.h"
#include <stdlib.h>

typedef struct node {
  int line_num;
  SXLINE *line;
  struct node *left;
  struct node *right;
} Node;

Node *new_node(int line_num, SXLINE *line) {
  Node *node = malloc(sizeof(Node));
  node->line = line;
  node->left = NULL;
  node->right = NULL;
  return node;
};

void insert_node(Node *root, int line_num, SXLINE *line) {
  if (root == NULL) {
    root = new_node(line_num, line);
    return;
  }
  if (root->line_num > line_num) {
    insert_node(root->left, line_num, line);
  } else {
    insert_node(root->right, line_num, line);
  }
}

Node* search_node(Node *root, int line_num) {
  if (root == NULL || root->line_num == line_num) {
    return root;
  }
  if (root->line_num > line_num) {
    return search_node(root->left, line_num);
  } else {
    return search_node(root->right, line_num);
  }
}

