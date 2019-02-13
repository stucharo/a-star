from calc import list_to_cls

s = list_to_cls()
while s.next is not None:
    s = s.next

print(s.x)
