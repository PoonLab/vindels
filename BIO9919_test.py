print('correct!')
def mascii(numlist):
    # map each number in the list `numlist` to a character the ASCII character set
    # i.e. convert each number to its corresponding ASCII character using the chr() function
    # and return the resulting string
    result = ''
    for i in numlist:
        result += chr(i)
    return result

input = [99, 111, 114, 114, 101, 99, 116, 33]
print(mascii(input))

d = {'A': 'c', 'B': 'o', 'C': 'r', 'D': 'e', 'E': 't', 'F': '!'}
input = 'feadccba'
output = ''

for i in range(len(input)-1, -1, -1):
    let = input[i]
    if let.islower():
        let = let.upper()
    output += d[let]

print(output)



letters = ['c', 'e', 'o', 'c', 'rt', 'r']
letters = "".join(letters)
unscrambled = ''
for offset in range(2):
    for i in range(offset, len(letters),2):
        unscrambled += letters[i]

print(unscrambled)
